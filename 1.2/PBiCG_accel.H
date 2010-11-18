/*
 * PBiCG_accel.H - part of the SpeedIT Classic toolkit
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

#ifndef PBiCG_ACCEL_H
#define PBiCG_ACCEL_H

/*
	Class allows to connect external BiCGstab solver written with CUDA to OpenFOAM. 
	
	NOTE: since many available solvers use Compressed Sparse Row matrix format,
	      class provides internal conversion from OpenFOAM matrix to CSR matrix.
*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PBiCG_accel Declaration
\*---------------------------------------------------------------------------*/

class PBiCG_accel
:
    public lduMatrix::solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        PBiCG_accel(const PBiCG_accel&);

        //- Disallow default bitwise assignment
        void operator=(const PBiCG_accel&);


public:

    //- Runtime type information
    TypeName("PBiCG_accel");


    // Constructors

        //- Construct from matrix components and solver data stream
        PBiCG_accel
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    // Destructor

        virtual ~PBiCG_accel()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual lduMatrix::solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




#endif 
