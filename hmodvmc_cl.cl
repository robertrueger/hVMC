/*
 * Copyright (c) 2012, Robert Rueger <rueger@itp.uni-frankfurt.de>
 *
 * This file is part of hVMC.
 *
 * hVMC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hVMC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hVMC.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef USE_FP_DBLPREC_OPENCL
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double fptype;
#else
typedef float  fptype;
#endif


kernel void update_W(
  global fptype* W_inbuf, global fptype* W_outbuf, uint W_cols,
  uint k, uint l, uint k_pos )
{
  // Note: storage order is row major
  // -> W_ij = W[ i * W_cols + j ]

  const uint i = get_global_id( 0 ) / W_cols;
  const uint j = get_global_id( 0 ) % W_cols;

  W_outbuf[ i * W_cols + j ] =
    fma(
      W_inbuf[     i * W_cols + k ] / W_inbuf[ l * W_cols + k ],
      W_inbuf[ k_pos * W_cols + j ] - W_inbuf[ l * W_cols + j ],
      W_inbuf[     i * W_cols + j ]
    );
}
