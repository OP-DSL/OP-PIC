/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include <opp_lib.h>
#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

#define FUSE_KERNELS

#define ONE                       1
#define TWO                       2
#define DIM                       2

#define VOXEL(x,y, nx0) (x + nx0 * y)

#define VOXEL_MAP(_ix,_iy, nx,ny, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    if (_x < 0) _x = (nx-1); \
    if (_y < 0) _y = (ny-1); \
    if (_x >= nx) _x = 0; \
    if (_y >= ny) _y = 0; \
    OUT = (_x + nx*_y); }

#define RANK_TO_INDEX(rank,ix,iy,_x) \
    int _ix, _iy;                                         \
    _ix  = (rank);        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(_x);   /* ix = ix */                   \
    (ix) = _ix;                                           \
    (iy) = _iy; 

enum Dim {
    x = 0,
    y = 1,
};

#define NEIGHBOURS 4

enum CellMap {  
    xd_y = 0,
    xu_y,
    x_yd,
    x_yu 
};

