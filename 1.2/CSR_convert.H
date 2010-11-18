/*
 * CSR_convert.H - part of the SpeedIT Classic toolkit
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
 
 #ifndef CSR_CONVERT_H
#define CSR_CONVERT_H

#include <vector>
#include "lduMatrix.H"


//
//	Conversion from OpenFOAM lduMatrix to classic Compressed Sparse Row matrix.
//
template<class REAL_TYPE>
void ldu2csr( 
							const Foam::lduMatrix           & matrix, 
										 std::vector< int >       & c_idx, 
										 std::vector< int >       & r_idx,
										 std::vector< REAL_TYPE > & vals
						) 
{
	int n_rows = matrix.diag().size() ;

	//
	//	Calculate each row size. Sizes are shifted, because in the next part
	//	array r_idx is modified twice.
	//
	for (int i=0; i < matrix.upper().size() ; i++) {
		int ri1 = matrix.lduAddr().lowerAddr()[i] ; 
		int ri2 = matrix.lduAddr().upperAddr()[i] ;
		r_idx[ri1+2] ++ ;
		r_idx[ri2+2] ++ ;
	} ;

	//
	//	Compute row offsets. Offsets are shifted by one positions, 
	//  because they are used as START positions while filling values.
	//
	r_idx[0] = 0 ;
	r_idx[1] = 0 ;
	for(int i=1; i<n_rows; i++) {
		r_idx[i+1] += r_idx[i] ;
	} ;

	//
	//	Fill in CSR matrix.
	//	Order below is important to keep column indices sorted.
	//

	// lower triangle
	for (int i=0; i < matrix.lower().size() ; i++) {
		int row    = matrix.lduAddr().upperAddr()[i] +1;
		int column = matrix.lduAddr().lowerAddr()[i] ;

		int idx = r_idx[row] ;
		vals[idx] = matrix.lower()[i] ;
		c_idx[idx] = column ;
		r_idx[row]++ ;
	} ;
	// diagonal
	for (int i=0; i<matrix.diag().size(); i++) {
		int idx = r_idx[i+1] ;
		vals[idx] = matrix.diag()[i] ;
		c_idx[idx] = i ; // i is row and column index
		r_idx[i+1]++ ;
	} ;
	// upper triangle
	for (int i=0; i < matrix.upper().size() ; i++) {
		int row    = matrix.lduAddr().lowerAddr()[i] +1;
		int column = matrix.lduAddr().upperAddr()[i] ;

		int idx = r_idx[row] ;
		vals[idx] = matrix.upper()[i] ;
		c_idx[idx] = column ;
		r_idx[row]++ ;
	} ;
} ;


#endif
