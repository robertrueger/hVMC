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


// ============================ PREPROCESSOR ===================================

#ifndef LATTICE_TYPE
# error "lattice type not defined"
#endif

#ifndef LATTICE_NUM_SITES
# error "number of lattice sites not defined"
#endif

#ifndef NUM_ELECTRONS
# error "number of electrons not defined"
#endif

#ifdef M_IS_SPIN_SYMMETRIC
# define W_ROWS LATTICE_NUM_SITES
# define W_COLS ( NUM_ELECTRONS / 2 )
#else
# define W_ROWS ( 2 * LATTICE_NUM_SITES )
# define W_COLS NUM_ELECTRONS
#endif



// ======================== DOUBLE PRECISION SETUP =============================

#ifdef USE_FP_DBLPREC_OPENCL
# pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double fptype;
#else
typedef float  fptype;
#endif



// ========================== LATTICE FUNCTIONS ================================


uint lattice_get_spinup_site( uint l )
{
  return l >= LATTICE_NUM_SITES ? l - LATTICE_NUM_SITES : l;
}

uint lattice_get_spinlinked_site( uint l )
{
  if ( l < LATTICE_NUM_SITES ) {
    return l + LATTICE_NUM_SITES;
  } else {
    return l - LATTICE_NUM_SITES;
  }
}


#if LATTICE_TYPE == 0 // =================== 1D CHAIN LATTICE ===================

constant uint LATTICE_NUM_NEIGHBORS[3] = { 2, 2, 2 };


uint lattice_get_1nn( uint l, uint id )
{
  if ( id == 0 ) { // left neighbor

    if ( l == 0 || l == LATTICE_NUM_SITES ) {
      return l + LATTICE_NUM_SITES - 1;
    } else {
      return l - 1;
    }

  } else { /* ( id == 1 ) */  // right neighbor

    if ( l == LATTICE_NUM_SITES - 1 || l == 2 * LATTICE_NUM_SITES - 1 ) {
      return l + 1 - LATTICE_NUM_SITES;
    } else {
      return l + 1;
    }

  }
}


uint lattice_get_2nn( uint l, uint id )
{
  if ( id == 0 ) { // 2nd left neighbor

    if ( l == 1 || l == LATTICE_NUM_SITES + 1 ) {
      return l + LATTICE_NUM_SITES - 2;
    } else {
      return l - 2;
    }

  } else { /* ( id == 1 ) */  // 2nd right neighbor

    if ( l == LATTICE_NUM_SITES - 2 || l == 2 * LATTICE_NUM_SITES - 2 ) {
      return l + 2 - LATTICE_NUM_SITES;
    } else {
      return l + 2;
    }

  }
}


uint lattice_get_3nn( uint l, uint id )
{
  if ( id == 0 ) { // 3rd left neighbor

    if ( l == 2 || l == LATTICE_NUM_SITES + 2 ) {
      return l + LATTICE_NUM_SITES - 3;
    } else {
      return l - 3;
    }

  } else { /* ( id == 1 ) */  // 3rd right neighbor

    if ( l == LATTICE_NUM_SITES - 3 || l == 2 * LATTICE_NUM_SITES - 3 ) {
      return l + 3 - LATTICE_NUM_SITES;
    } else {
      return l + 3;
    }

  }
}


uint lattice_reduce_idxrel( uint i, uint j )
{
  return min( abs_diff( i, j ), LATTICE_NUM_SITES );
}



#elif LATTICE_TYPE == 1 // ================ 2D SQUARE LATTICE ===================

#ifndef LATTICE_SIZE
#error "lattice size not defined for the 2d square lattice"
#endif

constant uint LATTICE_NUM_NEIGHBORS[3] = { 4, 4, 4 };


uint lattice_get_1nn( uint l, uint id )
{
  if ( id == 0 ) { // left neighbor

    if ( l % LATTICE_SIZE == 0 ) {
      return l + LATTICE_SIZE - 1;
    } else {
      return l - 1;
    }

  } else if ( id == 1 ) { // right neighbor

    if ( ( l + 1 ) % LATTICE_SIZE == 0 ) {
      return l + 1 - LATTICE_SIZE;
    } else {
      return l + 1;
    }

  } else if ( id == 2 ) { // bottom neighbor

    if ( l < LATTICE_SIZE ||
         ( l >= LATTICE_NUM_SITES && l < LATTICE_NUM_SITES + LATTICE_SIZE )
       ) {
      return l + LATTICE_NUM_SITES - LATTICE_SIZE;
    } else {
      return l - LATTICE_SIZE;
    }

  } else { /* ( id == 3 ) */  // top neighbor

    if ( l >= 2 * LATTICE_NUM_SITES - LATTICE_SIZE ||
         ( l < LATTICE_NUM_SITES && l >= LATTICE_NUM_SITES - LATTICE_SIZE )
       ) {
      return l + LATTICE_SIZE - LATTICE_NUM_SITES;
    } else {
      return l + LATTICE_SIZE;
    }

  }
}


uint lattice_get_2nn( uint l, uint id )
{
  if ( id == 0 ) { // bottom left neighbor

    if ( l % LATTICE_SIZE == 0 ) {
      // in left column
      if ( l == 0 || l == LATTICE_NUM_SITES ) {
        // bottom left corner
        return l + LATTICE_NUM_SITES - 1;
      } else {
        return l - 1;
      }
    } else if ( l < LATTICE_SIZE ||
                ( l >= LATTICE_NUM_SITES &&
                  l < LATTICE_NUM_SITES + LATTICE_SIZE
                )
              ) {
      // in bottom row
      // (but NOT in bottom left corner!)
      return l + LATTICE_NUM_SITES - LATTICE_SIZE - 1;
    } else {
      // in the center
      return l - LATTICE_SIZE - 1;
    }

  } else if ( id == 1 ) { // bottom right neighbor

    if ( ( l + 1 ) % LATTICE_SIZE == 0 ) {
      // in right column
      if ( l == LATTICE_SIZE - 1 || l == LATTICE_NUM_SITES + LATTICE_SIZE - 1 ) {
        // bottom right corner
        return l + LATTICE_NUM_SITES + 1 - 2 * LATTICE_SIZE;
      } else {
        return l + 1 - 2 * LATTICE_SIZE;
      }
    } else if ( l < LATTICE_SIZE ||
                ( l >= LATTICE_NUM_SITES &&
                  l < LATTICE_NUM_SITES + LATTICE_SIZE
                )
              ) {
      // in bottom row
      // (but NOT in bottom right corner!)
      return l + LATTICE_NUM_SITES + 1 - LATTICE_SIZE;
    } else {
      // in the center
      return l + 1 - LATTICE_SIZE;
    }

  } else if ( id == 2 ) { // top right neighbor

    if ( ( l + 1 ) % LATTICE_SIZE == 0 ) {
      // in right column
      if ( l == LATTICE_NUM_SITES - 1 || l == 2 * LATTICE_NUM_SITES - 1 ) {
        // top right corner
        return l + 1 - LATTICE_NUM_SITES;
      } else {
        return l + 1;
      }
    } else if ( l >= 2 * LATTICE_NUM_SITES - LATTICE_SIZE ||
                ( l < LATTICE_NUM_SITES &&
                  l >= LATTICE_NUM_SITES - LATTICE_SIZE
                )
              ) {
      // in top row
      // (but NOT in top right corner!)
      return l + LATTICE_SIZE + 1 - LATTICE_NUM_SITES;
    } else {
      // in the center
      return l + LATTICE_SIZE + 1;
    }

  } else { /* ( id == 3 ) */  // top left neighbor

    if ( l % LATTICE_SIZE == 0 ) {
      // in left column
      if ( l == LATTICE_NUM_SITES - LATTICE_SIZE ||
           l == 2 * LATTICE_NUM_SITES - LATTICE_SIZE
         ) {
        // top left corner
        return l + 2 * LATTICE_SIZE - 1 - LATTICE_NUM_SITES ;
      } else {
        return l + 2 * LATTICE_SIZE - 1;
      }
    } else if ( l >= 2 * LATTICE_NUM_SITES - LATTICE_SIZE ||
                ( l < LATTICE_NUM_SITES &&
                  l >= LATTICE_NUM_SITES - LATTICE_SIZE
                )
              ) {
      // in top row
      // (but NOT in top left corner!)
      return l + LATTICE_SIZE - 1 - LATTICE_NUM_SITES;
    } else {
      // in the center
      return l + LATTICE_SIZE - 1;
    }

  }
}


uint lattice_get_3nn( uint l, uint id )
{
  if ( id == 0 ) { // 2nd left neighbor

    if ( l % LATTICE_SIZE <= 1  ) {
      return l + LATTICE_SIZE - 2;
    } else {
      return l - 2;
    }

  } else if ( id == 1 ) { // 2nd right neighbor

    if ( l % LATTICE_SIZE >= LATTICE_SIZE - 2 ) {
      return l + 2 - LATTICE_SIZE;
    } else {
      return l + 2;
    }

  } else if ( id == 2 ) { // 2nd bottom neighbor

    if ( l < 2 * LATTICE_SIZE ||
         ( l >= LATTICE_NUM_SITES && l < LATTICE_NUM_SITES + 2 * LATTICE_SIZE )
       ) {
      return l + LATTICE_NUM_SITES - 2 * LATTICE_SIZE;
    } else {
      return l - 2 * LATTICE_SIZE;
    }

  } else { /* ( id == 3 ) */  // 2nd top neighbor

    if ( l >= 2 * ( LATTICE_NUM_SITES - LATTICE_SIZE ) ||
         ( l < LATTICE_NUM_SITES && l >= LATTICE_NUM_SITES - 2 * LATTICE_SIZE )
       ) {
      return l + 2 * LATTICE_SIZE - LATTICE_NUM_SITES;
    } else {
      return l + 2 * LATTICE_SIZE;
    }

  }
}


uint lattice_reduce_idxrel( uint i, uint j )
{
  // calculate the positions of i and j
  const uint x_i = i % LATTICE_SIZE;
  const uint y_i = i / LATTICE_SIZE;
  const uint x_j = j % LATTICE_SIZE;
  const uint y_j = j / LATTICE_SIZE;

  // calculate the position difference
  uint dx = abs_diff( x_i, x_j );
  uint dy = abs_diff( y_i, y_j );

  // wrap large differences around the boundaries
  if ( dx > LATTICE_SIZE / 2 ) {
    dx = LATTICE_SIZE - dx;
  }
  if ( dy > LATTICE_SIZE / 2 ) {
    dy = LATTICE_SIZE - dy;
  }

  // dx should be larger than dy
  if ( dy > dx ) {
    uint temp = dy;
    dy = dx;
    dx = temp;
  }

  return dx + LATTICE_SIZE * dy;
}

#undef LATTICE_SIZE // it is internal to the 2d square lattice

#else // LATTICE_TYPE != 0 or 1
#error "unknown lattice type"
#endif



// ======================= ELECTRON HOP DATA STRUCTURE =========================


struct ElectronHop {
  uint k; // id of the hopping electron
  uint l; // site that it hops to
  uint k_pos; // position of electron k before the hop
  bool accepted; // whether the hop has been accepted
};



// ============================== KERNELS ======================================

// Note: storage order is row major
// -> W_ij = W[ i * W_COLS + j ]


kernel void hop(
  global struct ElectronHop* hop_out,
  global uint* electron_pos, global uint* site_occ,
  global const fptype* Wbu, global const fptype* Wd,
  global const fptype* T, constant fptype* expv,
  uint rand_hopel_id, uint rand_hopnn_id, fptype rand_selector )
{
  if ( get_global_id( 0 ) == 0 ) {

#if VERBOSE
    printf( "Electron configuration is \n" );
    for ( uint i = 0; i < NUM_ELECTRONS; ++i ) {
      printf( "%i ", electron_pos[ i ] );
    }
    printf( "\n" );
    for ( uint l = 0; l < 2 * LATTICE_NUM_SITES; ++l ) {
      if ( l == LATTICE_NUM_SITES ) {
        printf( "\n" );
      }
      printf( "%i ", site_occ[ l ] );
    }
    printf( "\n" );
#endif

    // propose a hop
    struct ElectronHop thishop;
    thishop.k        = rand_hopel_id;
    thishop.k_pos    = electron_pos[ rand_hopel_id ];
    thishop.l        = lattice_get_1nn( thishop.k_pos, rand_hopnn_id );
    thishop.accepted = false;

#if VERBOSE
    printf(
      "Proposed hop: k = %i, l = %i, k_pos = %i\n",
      thishop.k, thishop.l, thishop.k_pos
    );
#endif

    // check if the hop is possible (= if site l is empty)
    if ( site_occ[ thishop.l ] == 0 ) {

#if VERBOSE
      printf( "Hop possible, site unoccupied!\n" );
#endif

      // calculate the acceptance probability
      const fptype R_j =
        T[ lattice_get_spinup_site( thishop.l ) ]
        / T[ lattice_get_spinup_site( thishop.k_pos ) ]
        * expv[ 0 ] / expv[ lattice_reduce_idxrel( thishop.l, thishop.k_pos ) ];

      const fptype R_s =
#ifdef M_IS_SPIN_SYMMETRIC
        ( thishop.k >= NUM_ELECTRONS / 2 ) ?
        Wd[ ( thishop.l - W_ROWS ) * W_COLS + ( thishop.k - W_COLS ) ] :
        Wbu[ thishop.l * W_COLS + thishop.k ];
#else
        Wbu[ thishop.l * W_COLS + thishop.k ];
#endif

      const fptype accept_prob = R_j * R_j * R_s * R_s;

      // check if the hop is accepted
      if ( accept_prob >= 1.f || rand_selector < accept_prob ) {

#if VERBOSE
        printf( "Hop accepted!\n" );
#endif

        // hop is accepted!
        thishop.accepted = true;

        // update the electronic configuration
        site_occ[ thishop.k_pos ] = 0;
        site_occ[ thishop.l ] = 1;
        electron_pos[ thishop.k ] = thishop.l;

#if VERBOSE
        printf( "New electron configuration is \n" );
        for ( uint i = 0; i < NUM_ELECTRONS; ++i ) {
          printf( "%i ", electron_pos[ i ] );
        }
        printf( "\n" );
        for ( uint l = 0; l < 2 * LATTICE_NUM_SITES; ++l ) {
          if ( l == LATTICE_NUM_SITES ) {
            printf( "\n" );
          }
          printf( "%i ", site_occ[ l ] );
        }
        printf( "\n" );
#endif

      }

#if VERBOSE
      else {
        printf( "Hop rejected!\n" );
      }
#endif

    }

#if VERBOSE
    else {
      printf( "Hop impossible, site is occupied!\n" );
    }
#endif

    // write the hop to global memory
    *hop_out = thishop;
  }
}



kernel void update_T(
  global struct ElectronHop* hop,
  global fptype* T_in, global fptype* T_out,
  constant fptype* expv )
{
  const uint i = get_global_id( 0 );

//  printf( "T_in[ %i ] = %f", i, T_in[ i ] );

  if ( hop->accepted ) {
    T_out[ i ] =
      T_in[ i ] *
      expv[ lattice_reduce_idxrel( i, lattice_get_spinup_site( hop->l ) ) ] /
      expv[ lattice_reduce_idxrel( i, lattice_get_spinup_site( hop->k_pos ) ) ];

//    printf( " -> T_out[ %i ] = %f\n", i, T_out[ i ] );
  } else {
    T_out[ i ] = T_in[ i ];
//    printf( " copied to T_out!\n" );
  }
}



kernel void update_Wbu(
  global const struct ElectronHop* hop,
  global const fptype* Wbu_in, global fptype* Wbu_out )
{
#ifdef M_IS_SPIN_SYMMETRIC
  // if there is spin symmetry in the single particle orbitals, this function
  // only needs to do something if the hopping electron has spin up!
  if ( hop->accepted && hop->l < LATTICE_NUM_SITES ) {
#else
  if ( hop->accepted ) {
#endif
    const uint i = get_global_id( 0 ) / W_COLS;
    const uint j = get_global_id( 0 ) % W_COLS;

    Wbu_out[ i * W_COLS + j ] =
      fma(

        Wbu_in[ i * W_COLS + hop->k ] / Wbu_in[ hop->l * W_COLS + hop->k ], // *

        Wbu_in[ hop->k_pos * W_COLS + j ] - Wbu_in[ hop->l * W_COLS + j ], // +

        Wbu_in[ i * W_COLS + j ]

      );
  } else {
    const uint k = get_global_id( 0 );
    Wbu_out[ k ] = Wbu_in[ k ];
  }
}



#ifdef M_IS_SPIN_SYMMETRIC
kernel void update_Wd(
  global const struct ElectronHop* hop,
  global const fptype* Wd_in, global fptype* Wd_out )
{
  if ( hop->accepted && hop->l >= LATTICE_NUM_SITES ) {
    const uint i = get_global_id( 0 ) / W_COLS;
    const uint j = get_global_id( 0 ) % W_COLS;

    Wd_out[ i * W_COLS + j ] =
      fma(

        Wd_in[ i * W_COLS + ( hop->k - W_COLS ) ]
        / Wd_in[ ( hop->l - W_ROWS ) * W_COLS + ( hop->k - W_COLS ) ], // *

        Wd_in[ ( hop->k_pos - W_ROWS ) * W_COLS + j ]
        - Wd_in[ ( hop->l - W_ROWS ) * W_COLS + j ], // +

        Wd_in[ i * W_COLS + j ]

      );
  } else {
    const uint k = get_global_id( 0 );
    Wd_out[ k ] = Wd_in[ k ];
  }
}
#endif



kernel void calc_E_l(
  global fptype* E_l_out,
  global const uint* electron_pos, global const uint* site_occ,
  global const fptype* Wbu, global const fptype* Wd,
  global const fptype* T, constant fptype* expv,
  constant fptype* Ut )
{
  const uint k = get_global_id( 0 );
  const uint k_pos = electron_pos[ k ];

  fptype E_l = 0.f;

  // 1: calculate kinetic energy part of E_l for electron k

  // loop over all neighbor orders X
  for ( uint X = 1; X <= 3; ++X ) {

    if ( Ut[ X ] == 0.f ) {
      continue;
    }

    fptype E_l_thisX = 0.f;

    // iterate over all neighbors Xni of order X
    for ( uint Xni = 0; Xni < LATTICE_NUM_NEIGHBORS[ X - 1 ]; ++Xni ) {

      // get the position l of the Xni'th neighbor of k_pos
      uint l;
      if ( X == 1 ) {
        l = lattice_get_1nn( k_pos, Xni );
      } else if ( X == 2 ) {
        l = lattice_get_2nn( k_pos, Xni );
      } else { /* X == 3 */
        l = lattice_get_3nn( k_pos, Xni );
      }

#if VERBOSE
      printf(
        "Considering hop of electron %i hopping from %i to %i (redidxrel is %i)\n",
        k, k_pos, l, lattice_reduce_idxrel( l, k_pos )
      );
#endif

      // check if electron k could hop to site l (it must be empty)
      if ( site_occ[ l ] == 0 ) {

        const fptype R_j = 1.f;
        T[ lattice_get_spinup_site( l ) ]
        / T[ lattice_get_spinup_site( k_pos ) ]
        * expv[ 0 ] / expv[ lattice_reduce_idxrel( l, k_pos ) ];

        const fptype R_s = 1.f;
#ifdef M_IS_SPIN_SYMMETRIC
        ( k >= NUM_ELECTRONS / 2 ) ?
        Wd[ ( l - W_ROWS ) * W_COLS + ( k - W_COLS ) ] :
        Wbu[ l * W_COLS + k ];
#else
        Wbu[ l * W_COLS + k ];
#endif

        E_l_thisX += R_j * R_s;
      }
    }

    E_l -= Ut[ X ] * E_l_thisX;
  }

  // 2: calculate potential energy part of E_l for electron k
  if ( k < NUM_ELECTRONS / 2 && site_occ[ k_pos + LATTICE_NUM_SITES ] == 1 ) {
    E_l += Ut[ 0 ];
  }

  // write results to the global memory
  E_l_out[ k ] = E_l;
}
