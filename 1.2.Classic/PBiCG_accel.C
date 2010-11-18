/*
 * PBiCG_accel.C - part of the SpeedIT Classic toolkit
 * Copyright 2010 (C) Vratis Ltd
 * email: support@vratis.com
 * 
 * SpeedIT Classic toolkit is a free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * SpeedIT Classic library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PBiCG_accel.H"
#include "CSR_convert.H"

#include "error.H"

//#include "si_classic.h"
#include "speedit_ex.h"

#include "cuda_runtime.h"

#include <vector>
#include <iostream>
#include <fstream>

// AK - MPI test
#include "mpi.h"
#include <cutil_inline.h>
#include <cutil.h>
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PBiCG_accel, 0);

    lduMatrix::solver::addasymMatrixConstructorToTable<PBiCG_accel>
        addPBiCG_accelAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBiCG_accel::PBiCG_accel
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//
//		Method "solve" is used as an interface to external solver functions
//
Foam::lduMatrix::solverPerformance Foam::PBiCG_accel::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{

	lduMatrix::solverPerformance solverPerf
		(
		 lduMatrix::preconditioner::getName(controlDict_) + typeName,
		 fieldName_
		);

	//
	// MPI multiGPU
	//
        
	int initialized;
	int myDevice, devMaxGflops;
	
	MPI_Initialized(&initialized);

	if(initialized)
	{
		int myrank, nprocs;
		int nDevices;
		int size_err, rank_err;

		devMaxGflops = cutGetMaxGflopsDeviceId();		
	
		size_err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        rank_err = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	

        cudaGetDeviceCount(&nDevices);
       
		int setDevice = (myrank+devMaxGflops) % nDevices; 
		cudaGetDevice(&myDevice);
		if(myDevice != setDevice) cudaSetDevice(setDevice);
	}
	else
	{
		devMaxGflops = cutGetMaxGflopsDeviceId();	
		cudaGetDevice(&myDevice);
		if(myDevice != devMaxGflops) cudaSetDevice(devMaxGflops);
	}
	
	//
	//	Create CSR matrix
	//
		
	int n_rows = matrix_.diag().size() ;
	int nnz = matrix_.lower().size() + matrix_.upper().size() + matrix_.diag().size() ;
	std::vector<scalar> vals(nnz) ;
	std::vector<int> c_idx(nnz) ;
	std::vector<int> r_idx(n_rows+2, 1) ;

	ldu2csr(matrix_, c_idx, r_idx, vals) ;

	scalar* p_vals = &(vals[0]) ;
	int* p_c_idx = &(c_idx[0]) ;
	int* p_r_idx = &(r_idx[0]) ;

	//
	// Get pointers for vectors containing solution and right hand side.
	//
	const scalar* B = source.cdata() ;
	scalar* X = psi.data() ;

	scalar eps = tolerance_ ;
	int n_iter = maxIter_ ;
	int solver_result = -1 ;

	//
	//	Choose preconditioner
	//
	PRECOND_TYPE precond = P_NONE ;
	
	if ("diagonal" == lduMatrix::preconditioner::getName(controlDict_)) {
		precond = P_DIAG ;
	} else {
		std::cout << "WARNING : Unsupported preconditioner " ;
		std::cout << lduMatrix::preconditioner::getName(controlDict_) ;
		std::cout << ", using NO preconditioner.\n" ;
	} ;

	//
	//	Copy data to GPU memory. Explicit memory magement is needed, because SpeedIT 
	//	Classic library solver functions require pointers to data in GPU memory.
	//
	scalar* pgpu_vals = NULL ;
	int* pgpu_c_idx   = NULL;
	int* pgpu_r_idx   = NULL ;
	scalar* pgpu_B    = NULL ;
	scalar* pgpu_X    = NULL ;

	int cu_err = 0 ;
	
	cu_err += cudaMalloc(&pgpu_vals,  nnz * sizeof(scalar)) ;
	cu_err += cudaMalloc(&pgpu_c_idx, nnz * sizeof(int)) ;
	cu_err += cudaMalloc(&pgpu_r_idx, (n_rows+1) * sizeof(int)) ;
	cu_err += cudaMalloc(&pgpu_X,     n_rows * sizeof(scalar)) ;
	cu_err += cudaMalloc(&pgpu_B,     n_rows * sizeof(scalar)) ;	

	if (0 != cu_err) {
		cudaFree (pgpu_vals) ;
		cudaFree (pgpu_c_idx) ;
		cudaFree (pgpu_B) ;
		cudaFree (pgpu_X) ;
		cudaFree (pgpu_r_idx) ;

		Foam::FatalError << "Can not allocate GPU memory" << Foam::exit(Foam::FatalError, -1) ;
	} ;

	cu_err += cudaMemcpy (pgpu_vals,  p_vals,  nnz        * sizeof(scalar), cudaMemcpyHostToDevice) ;
	cu_err += cudaMemcpy (pgpu_c_idx, p_c_idx, nnz        * sizeof(int),    cudaMemcpyHostToDevice) ;
	cu_err += cudaMemcpy (pgpu_B,     B,       n_rows     * sizeof(scalar), cudaMemcpyHostToDevice) ;
	cu_err += cudaMemcpy (pgpu_X,     X,       n_rows     * sizeof(scalar), cudaMemcpyHostToDevice) ;
	cu_err += cudaMemcpy (pgpu_r_idx, p_r_idx, (n_rows+1) * sizeof(int),    cudaMemcpyHostToDevice) ;

	if (0 != cu_err) {
		cudaFree (pgpu_vals) ;
		cudaFree (pgpu_c_idx) ;
		cudaFree (pgpu_B) ;
		cudaFree (pgpu_X) ;
		cudaFree (pgpu_r_idx) ;

		Foam::FatalError << "Can not copy data to GPU memory" << Foam::exit(Foam::FatalError, -1) ;
	} ;


		/*
							INSERT SOLVER FUNCTION HERE
		*/
#if defined(WM_SP)     // scalar is float

	Foam::FatalError << "Call your own solver from PBiCG_accel.C file" << Foam::exit(Foam::FatalError, -1) ;

	//solver_result = si_cscsrbicgstab( n_rows, 
	//								 									p_vals, p_c_idx, p_r_idx, 
	//											 			 			X, B,
	//																	precond,
	//								 			 						&n_iter, &eps) ;

	//solver_result = si_gscsrbicgstab( n_rows, 
	//							 									pgpu_vals, pgpu_c_idx, pgpu_r_idx, 
	//										 			 			pgpu_X, pgpu_B,
	//																precond,
	//							 			 						&n_iter, &eps) ;



#elif defined(WM_DP)   // scalar is double
		
	Foam::FatalError << "Call your own solver from PBiCG_accel.C file" << Foam::exit(Foam::FatalError, -1) ;

	//	solver_result = si_cdcsrbicgstab( n_rows, 
	//									 									p_vals, p_c_idx, p_r_idx, 
	//												 			 			X, B,
	//																		precond,
	//									 			 						&n_iter, &eps) ;

	//solver_result = si_gdcsrbicgstab( n_rows, 
	//								 								pgpu_vals, pgpu_c_idx, pgpu_r_idx, 
	//											 			 		pgpu_X, pgpu_B,
	//																precond,
	//								 			 					&n_iter, &eps) ;

#endif

	//
	//	Copy result from GPU memory
	//
	if (0 != cudaMemcpy (X, pgpu_X, n_rows*sizeof(scalar), cudaMemcpyDeviceToHost)) {
		std::cout << "ERROR : Can not copy result from GPU memory\n" ;
	} ;


	if (0 != solver_result) {
		//std::cout << "ERROR : " << si_errstr(solver_result) << "\n" ;
		std::cout << "ERROR : solver function returned " << solver_result << "\n" ;
	} ;

	//
	//	Free buffers in GPU memory
	//
	cudaFree (pgpu_vals) ;
	cudaFree (pgpu_c_idx) ;
	cudaFree (pgpu_B) ;
	cudaFree (pgpu_X) ;
	cudaFree (pgpu_r_idx) ;

	solverPerf.finalResidual() = eps ;
	solverPerf.nIterations() = n_iter ;
	solverPerf.checkConvergence(tolerance_, relTol_) ;

	return solverPerf ;
}


