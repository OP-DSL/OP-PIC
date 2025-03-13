
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

#include <hdf5.h>
#include <H5FDmpio.h>
#include "opp_lib.h"

opp_set opp_decl_set_hdf5(char const *file, char const *name);
opp_set opp_decl_particle_set_hdf5(char const *file, char const *name, opp_set cells_set);

opp_map opp_decl_map_hdf5(opp_set from, opp_set to, int dim, char const *file, char const *name);
opp_dat opp_decl_dat_hdf5(opp_set set, int dim, opp_data_type dtype, char const *file, char const *name);

void opp_get_const_hdf5(char const *name, int dim, char const *type, char *const_data, char const *file_name);

void opp_dump_to_hdf5(char const *file_name);

void opp_fetch_data_hdf5_file(opp_dat dat, char const *file_name);
void opp_fetch_data_hdf5_file_path(opp_dat dat, char const *file_name, char const *path_name);
