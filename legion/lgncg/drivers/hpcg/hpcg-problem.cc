/**
 * Copyright (c) 2014      Los Alamos National Security, LLC
 *                         All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * LA-CC 10-123
 */

#include "hpcg-problem.h"
#include "hpcg-problem-generator.h"
#include "hpcg-geometry.h"

int Problem::genProbTID       = -1;
int Problem::populatef2cTID   = -1;

// again - ideally, we just call into some central place to get TIDs.
// pretend this isn't here...
enum {
    GEN_PROB_TID = 2,
    POP_F2C_FID
};

namespace {
/**
 * task arguments for genProbTasks
 */
struct GPTaskArgs {
    // problem geometry
    Geometry geom;
    int64_t stencilSize;
    lgncg::SparseMatrix sa;
    lgncg::Vector va;
    lgncg::Vector vb;

    GPTaskArgs(void) : geom(0, 0, 0), stencilSize(0) { ; }

    GPTaskArgs(lgncg::SparseMatrix &sa,
               lgncg::Vector *va,
               lgncg::Vector *vb,
               Geometry geom,
               int64_t stencilSize)
    {
        // will always have a sparse matrix
        this->sa = sa;
        // only will have vectors sometimes
        if (va) this->va = *va;
        if (vb) this->vb = *vb;
        this->geom = geom;
        this->stencilSize = stencilSize;
    }
};

/**
 * task arguments for populatef2cTasks
 */
struct PF2CTaskArgs {
    // fine to coarse operator vector
    lgncg::Vector f2c;
    // fine geometry
    Geometry fGeom;
    // coarse geometry
    Geometry cGeom;

    PF2CTaskArgs(const lgncg::Vector &f2c,
                 const Geometry &fGeom,
                 const Geometry cGeom)
        : f2c(f2c), fGeom(fGeom), cGeom(cGeom) { ; }
};

} // end namespace

/**
 * wrapper to launch genProblem. always fill from left to right. A -> x -> y.
 */
void
Problem::genProb(lgncg::SparseMatrix &A,
                 lgncg::Vector *x,
                 lgncg::Vector *b,
                 const Geometry &geom,
                 int64_t stencilSize,
                 LegionRuntime::HighLevel::Context &ctx,
                 LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using LegionRuntime::Arrays::Rect;

    int idx = 0;
    // setup task args
    GPTaskArgs tArgs(A, x, b, geom, stencilSize);
    // setup task launcher
    TaskLauncher tl(Problem::genProbTID, TaskArgument(&tArgs, sizeof(tArgs)));
    // for each logical region
    //     add a region requirement
    //     and for each field the region contains
    //         add it to the launcher
    // WRITE_DISCARD is a special form of READ_WRITE that still permits
    // the task to perform any kind of operation, but informs the
    // runtime that the task intends to overwrite all previous data
    // stored in the logical region without reading it.  init launcher
    // A's vals
    tl.add_region_requirement(
        RegionRequirement(A.vals.lr, WRITE_DISCARD, EXCLUSIVE, A.vals.lr)
    );
    tl.add_field(idx++, A.vals.fid);
    // A's diag values
    tl.add_region_requirement(
        RegionRequirement(A.diag.lr, WRITE_DISCARD, EXCLUSIVE, A.diag.lr)
    );
    tl.add_field(idx++, A.diag.fid);
    // A's mIdxs
    tl.add_region_requirement(
        RegionRequirement(A.mIdxs.lr, WRITE_DISCARD, EXCLUSIVE, A.mIdxs.lr)
    );
    tl.add_field(idx++, A.mIdxs.fid);
    // A's # non 0s in row
    tl.add_region_requirement(
        RegionRequirement(A.nzir.lr, WRITE_DISCARD, EXCLUSIVE, A.nzir.lr)
    );
    tl.add_field(idx++, A.nzir.fid);
    // if we are provided a pointer to an x, add it
    if (x) {
        tl.add_region_requirement(
            RegionRequirement(x->lr, WRITE_DISCARD, EXCLUSIVE, x->lr)
        );
        tl.add_field(idx++, x->fid);
    }
    // similarly for b
    if (b) {
        tl.add_region_requirement(
            RegionRequirement(b->lr, WRITE_DISCARD, EXCLUSIVE, b->lr)
        );
        tl.add_field(idx++, b->fid);
    }
    // ... and go!
    (void)lrt->execute_task(ctx, tl);
}

/**
 * akin to HPCG's GenerateProblem.
 */
void
genProbTask(
    const LegionRuntime::HighLevel::Task *task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using namespace lgncg;
    // first capture how many regions were passed to us. this number dictates
    // what regions we will be working on. region order here matters.
    size_t nRgns = rgns.size();
    // we always have 4 regions for the sparse matrix for this task
    assert(nRgns >= 4 && nRgns <= 6);
    bool haveX = nRgns > 4;
    bool haveB = nRgns > 5;
    size_t rid = 0;
    GPTaskArgs targs = *(GPTaskArgs *)task->args;
    ////////////////////////////////////////////////////////////////////////////
    // Sparse Matrix A
    ////////////////////////////////////////////////////////////////////////////
    // remember that sparse matrix A has three regions
    const PhysicalRegion &avpr = rgns[rid++];
    const PhysicalRegion &adpr = rgns[rid++];
    const PhysicalRegion &aipr = rgns[rid++];
    const PhysicalRegion &azpr = rgns[rid++];
    // construct new initial conditions
    ProblemGenerator ic(targs.geom.nx,
                        targs.geom.ny,
                        targs.geom.nz,
                        targs.stencilSize);
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>  GDRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t> GLRA;
    typedef RegionAccessor<AccessorType::Generic, uint8_t> GSRA;
    // get handles to all the matrix accessors that we need
    GDRA aVals = avpr.get_field_accessor(targs.sa.vals.fid).typeify<double>();
    GDRA aDiag = adpr.get_field_accessor(targs.sa.diag.fid).typeify<double>();
    GLRA aMidx = aipr.get_field_accessor(targs.sa.mIdxs.fid).typeify<int64_t>();
    GSRA aNZir = azpr.get_field_accessor(targs.sa.nzir.fid).typeify<uint8_t>();

    IndexSpace indxspc = avpr.get_logical_region().get_index_space();
    IndexSpace dindxspc = adpr.get_logical_region().get_index_space();
    // the bounds of vals and midx
    Rect<1> pbnds = lrt->get_index_space_domain(ctx, indxspc).get_rect<1>();
    Rect<1> qbnds = lrt->get_index_space_domain(ctx, dindxspc).get_rect<1>();
    typedef GenericPointInRectIterator<1> GPRI1D;
    typedef DomainPoint DomPt;
    { // aVals and aMidx
        GPRI1D p(pbnds);
        GPRI1D q(qbnds);
        for (int64_t i = 0; i < targs.sa.nRows; ++i, q++) {
            for (int64_t j = 0; j < targs.sa.nCols; ++j, p++) {
                aVals.write(DomPt::from_point<1>(p.p), ic.A[i][j]);
                aMidx.write(DomPt::from_point<1>(p.p), ic.mtxInd[i][j]);
            }
            // SKG - i'm a little confused by why matDiag is not just an array
            aDiag.write(DomPt::from_point<1>(q.p), ic.matDiag[i][0]);
        }
    }
    indxspc = azpr.get_logical_region().get_index_space();
    pbnds = lrt->get_index_space_domain(ctx, indxspc).get_rect<1>();
    { // nzir
        GPRI1D p(pbnds);
        for (int64_t i = 0; p; ++i, p++) {
            aNZir.write(DomPt::from_point<1>(p.p), ic.non0sInRow[i]);
        }
    }
    if (haveX) {
        ////////////////////////////////////////////////////////////////////////
        // Vector x - initial guess of all zeros
        ////////////////////////////////////////////////////////////////////////
        const PhysicalRegion &vecXPr  = rgns[rid++];
        GDRA vVals = vecXPr.get_field_accessor(targs.va.fid).typeify<double>();
        indxspc = vecXPr.get_logical_region().get_index_space();
        pbnds = lrt->get_index_space_domain(ctx, indxspc).get_rect<1>();
        for (GPRI1D p(pbnds); p; p++) {
            vVals.write(DomPt::from_point<1>(p.p), 0.0);
        }
    }
    if (haveB) {
        ////////////////////////////////////////////////////////////////////////
        // Vector b
        ////////////////////////////////////////////////////////////////////////
        const PhysicalRegion &vecBPr  = rgns[rid++];
        GDRA vVals = vecBPr.get_field_accessor(targs.vb.fid).typeify<double>();
        indxspc = vecBPr.get_logical_region().get_index_space();
        pbnds = lrt->get_index_space_domain(ctx, indxspc).get_rect<1>();
        GPRI1D p(pbnds);
        for (int64_t i = 0; p; ++i, p++) {
            vVals.write(DomPt::from_point<1>(p.p), ic.b[i]);
        }
    }
}

/**
 * wrapper to launch genCoarseProblem.
 */
void
Problem::genCoarseProb(lgncg::SparseMatrix &Af,
                       const Geometry &fineGeom,
                       Geometry &coarseGeom,
                       int64_t stencilSize,
                       LegionRuntime::HighLevel::Context &ctx,
                       LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using LegionRuntime::Arrays::Rect;

    int64_t nx = fineGeom.nx;
    int64_t ny = fineGeom.ny;
    int64_t nz = fineGeom.nz;
    // sanity
    assert(0 == (nx % 2) && 0 == (ny % 2) && 0 == (nz % 2));
    assert(Af.nRows == nx * ny * nz);
    // set the coarse geometry
    coarseGeom = Geometry(nx / 2, ny / 2, nz / 2);
    // now construct the required data structures for the coarse grid
    Af.Ac = new lgncg::SparseMatrix();
    assert(Af.Ac);
    int64_t cnx = nx / 2;
    int64_t cny = ny / 2;
    int64_t cnz = nz / 2;
    // now create the coarse grid matrix
    Af.Ac->create(cnx * cny * cnz, Af.nCols, ctx, lrt);
    Af.mgData = new lgncg::MGData(nx  * ny  * nz,
                                  cnx * cny * cnz,
                                  ctx, lrt);
    assert(Af.mgData);
    // start task launch prep
    int idx = 0;
    // setup task launcher to populate f2cOp
    PF2CTaskArgs targs(Af.mgData->f2cOp, fineGeom, coarseGeom);
    TaskLauncher tl(Problem::populatef2cTID,
                    TaskArgument(&(targs), sizeof(targs)));
    lgncg::MGData *mgd = Af.mgData;
    // f2cOp
    tl.add_region_requirement(
        RegionRequirement(mgd->f2cOp.lr, WRITE_DISCARD, EXCLUSIVE, mgd->f2cOp.lr)
    );
    tl.add_field(idx++, mgd->f2cOp.fid);
    // populate the fine to coarse operator vector
    lrt->execute_task(ctx, tl).get_void_result();
    // now generate problem at this level
    genProb(*Af.Ac, NULL, NULL, coarseGeom, stencilSize, ctx, lrt);
}

/**
 * this task is responsible for populating the fine to coarse operator
 */
void
populatef2cTask(
    const LegionRuntime::HighLevel::Task *task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using namespace lgncg;
    // f2cOp
    assert(1 == rgns.size());
    size_t rid = 0;
    PF2CTaskArgs args = *(PF2CTaskArgs *)task->args;
    const lgncg::Vector &f2cv = args.f2c;
    const Geometry &fGeom = args.fGeom;
    const Geometry &cGeom = args.cGeom;
    // name the regions
    const PhysicalRegion &f2cpr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, int64_t>  GIRA;
    // get handles to all the matrix accessors that we need
    // fine length (both are the same)
    GIRA f2c = f2cpr.get_field_accessor(f2cv.fid).typeify<int64_t>();
    // setup the index spaces that we need.
    IndexSpace findxspc = f2cpr.get_logical_region().get_index_space();
    // similarly, do the same for the bounds
    Rect<1> fbnds = lrt->get_index_space_domain(ctx, findxspc).get_rect<1>();

    // TODO - add offset density checks -- maybe?
    Rect<1> f2csr; ByteOffset f2cOff[1];
    int64_t *f2cp = f2c.raw_rect_ptr<1>(f2cv.bounds, f2csr, f2cOff);

    for (int64_t i = 0; i < fbnds.volume(); ++i) {
        f2cp[i] = 0;
    }
    int64_t nxf = fGeom.nx;
    int64_t nyf = fGeom.ny;
    int64_t nxc = cGeom.nx;
    int64_t nyc = cGeom.ny;
    int64_t nzc = cGeom.nz;
    for (int64_t izc = 0; izc < nzc; ++izc) {
        int64_t izf = 2 * izc;
        for (int64_t iyc = 0; iyc < nyc; ++iyc) {
            int64_t iyf = 2 * iyc;
            for (int64_t ixc = 0; ixc < nxc; ++ixc) {
                int64_t ixf = 2 * ixc;
                int64_t currentCoarseRow = izc * nxc * nyc + iyc * nxc + ixc;
                int64_t currentFineRow = izf * nxf * nyf + iyf * nxf + ixf;
                f2cp[currentCoarseRow] = currentFineRow;
            } // end iy loop
        } // end even iz if statement
    } // end iz loop
}

/**
 * registers any legion tasks that Problem may use.
 */
void
Problem::init(void)
{
    using namespace LegionRuntime::HighLevel;
    // register any tasks here -- exactly once
    static bool registeredTasks = false;
    if (!registeredTasks) {
        Problem::genProbTID = GEN_PROB_TID;
        // register the generate problem task
        HighLevelRuntime::register_legion_task<genProbTask>(
            Problem::genProbTID /* task id */,
            Processor::LOC_PROC /* proc kind  */,
            true /* single */,
            false /* index */,
            AUTO_GENERATE_ID,
            TaskConfigOptions(false /* leaf task */),
            "gen-prob-task"
        );
        // register the populate f2cOperator task
        Problem::populatef2cTID = POP_F2C_FID;
        HighLevelRuntime::register_legion_task<populatef2cTask>(
            Problem::populatef2cTID /* task id */,
            Processor::LOC_PROC /* proc kind  */,
            true /* single */,
            false /* index */,
            AUTO_GENERATE_ID,
            TaskConfigOptions(false /* leaf task */),
            "populate-f2c-task"
        );
        registeredTasks = true;
    }
}

