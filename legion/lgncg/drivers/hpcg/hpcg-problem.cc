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

#include "lgncg.h"

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
 * task arguments for HPCG driver
 */
struct HPCGTaskArgs {
    lgncg::DSparseMatrix sa;
    lgncg::DVector va;
    lgncg::DVector vb;

    HPCGTaskArgs(const lgncg::SparseMatrix &sma,
                 const lgncg::Vector *const vpa,
                 const lgncg::Vector *const vpb) {
        sa = sma;
        if (vpa) va = *vpa;
        if (vpb) vb = *vpb;
    }
};

/**
 * task arguments for populatef2cTasks
 */
struct PF2CTaskArgs {
    // fine to coarse operator vector
    lgncg::Vector f2c;
    // fine geometry
    lgncg::Geometry fGeom;
    // coarse geometry
    lgncg::Geometry cGeom;

    PF2CTaskArgs(const lgncg::Vector &f2c,
                 const lgncg::Geometry &fGeom,
                 const lgncg::Geometry cGeom)
        : f2c(f2c), fGeom(fGeom), cGeom(cGeom) { ; }
};

} // end namespace

/**
 * wrapper to launch setICsTask. always fill from left to right. A -> x -> y.
 */
void
Problem::setICs(lgncg::SparseMatrix &A,
                lgncg::Vector *x,
                lgncg::Vector *b,
                LegionRuntime::HighLevel::Context &ctx,
                LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using LegionRuntime::Arrays::Rect;

    int idx = 0;
    // setup task args
    HPCGTaskArgs targs(A, x, b);
    ArgumentMap argMap;
    const int64_t npx = A.geom.npx;
    const int64_t npy = A.geom.npy;
    for (int i = 0; i < A.vals.lDom().get_volume(); ++i) {
        targs.sa.vals.sgb  = A.vals.sgb()[i];
        targs.sa.diag.sgb  = A.diag.sgb()[i];
        targs.sa.mIdxs.sgb = A.mIdxs.sgb()[i];
        targs.sa.nzir.sgb  = A.nzir.sgb()[i];
        targs.sa.g2g.sgb   = A.g2g.sgb()[i];
        // setup task to global geometry
        targs.sa.geom.ipz = i / (npx * npy);
        targs.sa.geom.ipy = (i - targs.sa.geom.ipz * npx * npy) / npx;
        targs.sa.geom.ipx = i % npx;
#if 0 // debug
        std::cout << "task: " << i
                  << " ipx: " << targs.sa.geom.ipx
                  << " ipy: " << targs.sa.geom.ipy
                  << " ipz: " << targs.sa.geom.ipz << std::endl;
#endif
        if (x) targs.va.sgb = x->sgb()[i];
        if (b) targs.vb.sgb = b->sgb()[i];
        argMap.set_point(DomainPoint::from_point<1>(Point<1>(i)),
                         TaskArgument(&targs, sizeof(targs)));
    }
    // setup task launcher
    IndexLauncher il(Problem::genProbTID, A.vals.lDom(),
                     TaskArgument(NULL, 0), argMap);
    // for each logical region
    //     add a region requirement
    //     and for each field the region contains
    //         add it to the launcher
    // WRITE_DISCARD is a special form of READ_WRITE that still permits
    // the task to perform any kind of operation, but informs the
    // runtime that the task intends to overwrite all previous data
    // stored in the logical region without reading it.  init launcher
    // A's vals
    il.add_region_requirement(
        RegionRequirement(A.vals.lp(), 0, WRITE_DISCARD, EXCLUSIVE, A.vals.lr)
    );
    il.add_field(idx++, A.vals.fid);
    // A's diag values
    il.add_region_requirement(
        RegionRequirement(A.diag.lp(), 0, WRITE_DISCARD, EXCLUSIVE, A.diag.lr)
    );
    il.add_field(idx++, A.diag.fid);
    // A's mIdxs
    il.add_region_requirement(
        RegionRequirement(A.mIdxs.lp(), 0, WRITE_DISCARD, EXCLUSIVE, A.mIdxs.lr)
    );
    il.add_field(idx++, A.mIdxs.fid);
    // A's # non 0s in row
    il.add_region_requirement(
        RegionRequirement(A.nzir.lp(), 0, WRITE_DISCARD, EXCLUSIVE, A.nzir.lr)
    );
    il.add_field(idx++, A.nzir.fid);
    // A's global to global table 
    il.add_region_requirement(
        RegionRequirement(A.g2g.lp(), 0, WRITE_DISCARD, EXCLUSIVE, A.g2g.lr)
    );
    il.add_field(idx++, A.g2g.fid);
    // if we are provided a pointer to an x, add it
    if (x) {
        il.add_region_requirement(
            RegionRequirement(x->lp(), 0, WRITE_DISCARD, EXCLUSIVE, x->lr)
        );
        il.add_field(idx++, x->fid);
    }
    // similarly for b
    if (b) {
        il.add_region_requirement(
            RegionRequirement(b->lp(), 0, WRITE_DISCARD, EXCLUSIVE, b->lr)
        );
        il.add_field(idx++, b->fid);
    }
    // ... and go!
    (void)lrt->execute_index_space(ctx, il);
}

/**
 * akin to HPCG's GenerateProblem.
 */
void
setICsTask(
    const LegionRuntime::HighLevel::Task *task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &rgns,
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *lrt
)
{
    using namespace LegionRuntime::HighLevel;
    using namespace LegionRuntime::Accessor;
    using namespace lgncg;

    (void)ctx;
    (void)lrt;
    // stash my task ID
    const int taskID = task->index_point.point_data[0];
    // capture how many regions were passed to us. this number dictates what
    // regions we will be working on. region order here matters.
    const size_t nRgns = rgns.size();
    // we always have 5 regions for the sparse matrix for this task
    assert(nRgns >= 5 && nRgns <= 7);
    const bool haveX = nRgns > 5;
    const bool haveB = nRgns > 6;
    size_t rid = 0;
    // get our task arguments
    const HPCGTaskArgs targs = *(HPCGTaskArgs *)task->local_args;
    ////////////////////////////////////////////////////////////////////////////
    // Sparse Matrix A
    ////////////////////////////////////////////////////////////////////////////
    const PhysicalRegion &avpr = rgns[rid++];
    const PhysicalRegion &adpr = rgns[rid++];
    const PhysicalRegion &aipr = rgns[rid++];
    const PhysicalRegion &azpr = rgns[rid++];
    const PhysicalRegion &g2gpr = rgns[rid++];
    // convenience typedefs
    typedef RegionAccessor<AccessorType::Generic, double>   GDRA;
    typedef RegionAccessor<AccessorType::Generic, int64_t>  GLRA;
    typedef RegionAccessor<AccessorType::Generic, uint8_t>  GSRA;
    typedef RegionAccessor<AccessorType::Generic, I64Tuple> GTRA;
    // get handles to all the matrix accessors that we need
    GDRA av = avpr.get_field_accessor(targs.sa.vals.fid).typeify<double>();
    GDRA ad = adpr.get_field_accessor(targs.sa.diag.fid).typeify<double>();
    GLRA ai = aipr.get_field_accessor(targs.sa.mIdxs.fid).typeify<int64_t>();
    GSRA az = azpr.get_field_accessor(targs.sa.nzir.fid).typeify<uint8_t>();
    GTRA at = g2gpr.get_field_accessor(targs.sa.g2g.fid).typeify<I64Tuple>();
    ////////////////////////////////////////////////////////////////////////////
    // all problem setup logic in the ProblemGenerator /////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // construct new initial conditions for this sub-region
    ProblemGenerator ic(targs.sa, taskID);
    // okay - now write out to regions. FIXME: just pass raw pointers to problem
    // setup to avoid memory bloat during init
    // the bounds of all entries
    typedef GenericPointInRectIterator<1> GPRI1D;
    typedef DomainPoint DomPt;
    do {
        const int64_t nLocalCols = targs.sa.nCols;
        const int64_t nLocalRows = targs.sa.diag.sgb.volume();
        GPRI1D p(targs.sa.vals.sgb);
        GPRI1D q(targs.sa.diag.sgb);
        for (int64_t i = 0; i < nLocalRows; ++i, q++) {
            for (int64_t j = 0; j < nLocalCols; ++j, p++) {
                av.write(DomPt::from_point<1>(p.p), ic.A[i][j]);
                ai.write(DomPt::from_point<1>(p.p), ic.mtxInd[i][j]);
            }
            // SKG - i'm a little confused by why matDiag is not just an array
            ad.write(DomPt::from_point<1>(q.p), ic.matDiag[i][0]);
            az.write(DomPt::from_point<1>(q.p), ic.non0sInRow[i]);
        }
    } while(0);
    // Populate g2g
    do {
        const int64_t offset = taskID * targs.sa.g2g.sgb.volume();
        int64_t row = 0;
        for (GPRI1D p(targs.sa.g2g.sgb); p; p++, row++) {
            at.write(DomPt::from_point<1>(p.p),
                     I64Tuple(offset + row, ic.l2gTab[row]));
        }
    } while(0);
    if (haveX) {
        ////////////////////////////////////////////////////////////////////////
        // Vector x - initial guess of all zeros
        ////////////////////////////////////////////////////////////////////////
        const PhysicalRegion &vecXPr = rgns[rid++];
        GDRA vVals = vecXPr.get_field_accessor(targs.va.fid).typeify<double>();
        for (GPRI1D p(targs.va.sgb); p; p++) {
            vVals.write(DomPt::from_point<1>(p.p), 0.0);
        }
    }
    if (haveB) {
        ////////////////////////////////////////////////////////////////////////
        // Vector b
        ////////////////////////////////////////////////////////////////////////
        const PhysicalRegion &vecBPr = rgns[rid++];
        GDRA vVals = vecBPr.get_field_accessor(targs.vb.fid).typeify<double>();
        GPRI1D p(targs.vb.sgb);
        for (int64_t i = 0; p; ++i, p++) {
            vVals.write(DomPt::from_point<1>(p.p), ic.b[i]);
        }
    }
}

/**
 * wrapper to launch genCoarseProblem.
 */
void
Problem::genCoarseProbGeom(lgncg::SparseMatrix &Af,
                           LegionRuntime::HighLevel::Context &ctx,
                           LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using namespace lgncg;
    using lgncg::Geometry;

    const int64_t npx = Af.geom.npx;
    const int64_t npy = Af.geom.npy;
    const int64_t npz = Af.geom.npz;
    const int64_t nx  = Af.geom.nx;
    const int64_t ny  = Af.geom.ny;
    const int64_t nz  = Af.geom.nz;
    const int64_t cnx = nx / 2;
    const int64_t cny = ny / 2;
    const int64_t cnz = nz / 2;
    const int64_t globalXYZ =  (npx * nx)  * (npy * ny)  * (npz * nz);
    const int64_t cGlobalXYZ = (npx * cnx) * (npy * cny) * (npz * cnz);
    // sanity
    assert(0 == (nx % 2) && 0 == (ny % 2) && 0 == (nz % 2));
    assert(Af.nRows == globalXYZ); 
    // now construct the required data structures for the coarse grid
    Af.Ac = new lgncg::SparseMatrix();
    assert(Af.Ac);
    // set the coarse geometry
    const Geometry coarseGeom = Geometry(Af.geom.size,
                                         npx, npy, npz,
                                         cnx, cny, cnz);
    // now create the coarse grid matrix
    Af.Ac->create(coarseGeom, cGlobalXYZ, Af.nCols, ctx, lrt);
    Af.Ac->partition(Af.geom, ctx, lrt);
    Af.mgData = new lgncg::MGData(globalXYZ, cGlobalXYZ, ctx, lrt);
    assert(Af.mgData);
    Af.mgData->partition(Af.nParts, ctx, lrt);
#if 0
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
    lrt->execute_task(ctx, tl);
#endif
}

#if 0
/**
 * wrapper to launch populatef2cTasks
 */
void
Problem::populatef2c(lgncg::SparseMatrix &Af,
                     const Geometry &fineGeom,
                     Geometry &coarseGeom,
                     LegionRuntime::HighLevel::Context &ctx,
                     LegionRuntime::HighLevel::HighLevelRuntime *lrt)
{
    using namespace LegionRuntime::HighLevel;
    using LegionRuntime::Arrays::Rect;
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
    (void)lrt->execute_task(ctx, tl);
#if 0
    // now generate problem at this level
    genProb(*Af.Ac, NULL, NULL, coarseGeom, stencilSize, ctx, lrt);
#endif
}
#endif

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
    const int64_t nxf = fGeom.npx * fGeom.nx;
    const int64_t nyf = fGeom.npy * fGeom.ny;
    const int64_t nxc = cGeom.npx * cGeom.nx;
    const int64_t nyc = cGeom.npy * cGeom.ny;
    const int64_t nzc = cGeom.npz * cGeom.nz;
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
        HighLevelRuntime::register_legion_task<setICsTask>(
            Problem::genProbTID /* task id */,
            Processor::LOC_PROC /* proc kind  */,
            false /* single */,
            true /* index */,
            AUTO_GENERATE_ID,
            TaskConfigOptions(true /* leaf task */),
            "set-ics-task"
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

