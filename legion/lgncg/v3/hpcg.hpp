
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file hpcg.hpp

 HPCG data structures and functions
 */

#ifndef HPCG_HPP
#define HPCG_HPP

#include "LegionStuff.hpp"

extern LegionRuntime::Logger::Category Logger;
#define HPCG_fout Logger.print()

struct HPCG_Params_STRUCT {
  int comm_size; //!< Number of MPI processes in MPI_COMM_WORLD
#if 0
  int comm_rank; //!< This process' MPI rank in the range [0 to comm_size - 1]
#endif
  int numThreads; //!< This process' number of threads
  int nx; //!< Number of x-direction grid points for each local subdomain
  int ny; //!< Number of y-direction grid points for each local subdomain
  int nz; //!< Number of z-direction grid points for each local subdomain
  int runningTime; //!< Number of seconds to run the timed portion of the benchmark
  int stencilSize; //!< Size of the stencil
};
/*!
  HPCG_Params is a shorthand for HPCG_Params_STRUCT
 */
typedef HPCG_Params_STRUCT HPCG_Params;

extern int HPCG_Init(HPCG_Params &params, const SPMDMeta &spmdMeta);
extern int HPCG_Finalize(void);

#endif // HPCG_HPP
