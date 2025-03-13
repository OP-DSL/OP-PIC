

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

#include <opp_mpi.h>

/*******************************************************************************
 * Main MPI halo creation routine
 *******************************************************************************/
void opp_halo_create() 
{
    __opp_halo_create();
}

/*******************************************************************************
 * Routine to Clean-up all MPI halos(called at the end of an OP2 MPI application)
*******************************************************************************/
void opp_halo_destroy() 
{
    __opp_halo_destroy();
}

/*******************************************************************************
 * Routine to start exchanging halos for the args
*******************************************************************************/
int opp_mpi_halo_exchanges_grouped(opp_set set, int nargs, opp_arg *args, DeviceType device)
{
    return opp_mpi_halo_exchanges(set, nargs, args); // Halo exchange only for mesh dats
}

/*******************************************************************************
 * Routine to exchange MPI halos of the all the args
*******************************************************************************/
int opp_mpi_halo_exchanges(opp_set set, int nargs, opp_arg *args) 
{
    if (OPP_DBG) opp_printf("opp_mpi_halo_exchanges", "START");

    int size = set->size;
    bool direct_flag = true;

    // check if this is a direct loop
    for (int n = 0; n < nargs; n++)
        if (args[n].opt && args[n].argtype == OPP_ARG_DAT && args[n].idx != -1)
            direct_flag = false;

    // return set size if it is a direct loop
    if (direct_flag)
        return size;

    // not a direct loop ...
    int exec_flag = 0;
    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].idx != -1 && args[n].acc != OPP_READ) 
        {
            size = set->size + set->exec_size;
            exec_flag = 1;
        }
    }

    for (int n = 0; n < nargs; n++) 
    {
        if (args[n].opt && args[n].argtype == OPP_ARG_DAT && (!args[n].dat->set->is_particle) && 
            (args[n].dat->dirtybit == 1)) 
        {
            bool already_done = false;

            // Check if dat was already done within these args
            for (int m = 0; m < n; m++) 
            {
                if ((args[n].dat == args[m].dat) && (args[m].idx != -1))
                    already_done = true;
            }

            if (!already_done)
            {
                if (OPP_DBG) 
                    opp_printf("opp_mpi_halo_exchanges", "opp_exchange_halo for dat [%s] exec_flag %d idx=%d", 
                        args[n].dat->name, exec_flag, args[n].idx);
                opp_mpi_halo_exchange(&args[n], exec_flag);
            }
        }
    }  

    return size;
}

/*******************************************************************************
 * Routine to exchange MPI halos of the arg
*******************************************************************************/
void opp_mpi_halo_exchange(opp_arg *arg, int exec_flag)
{
    if (arg->opt == 0)
        return;

    // For a directly accessed opp_dat do not do halo exchanges if not executing over redundant compute block
    if (exec_flag == 0 && arg->idx == -1)
        return;

    opp_dat dat = arg->dat;

        // need to exchange both direct and indirect data sets if they are dirty
    if ((arg->acc == OPP_READ || arg->acc == OPP_RW) && (dat->dirtybit == 1)) 
    {
        if (OPP_DBG) 
            opp_printf("opp_mpi_halo_exchange", "__opp_mpi_host_halo_exchange for dat [%s] exec_flag %d", 
                arg->dat->name, exec_flag);
        __opp_mpi_host_halo_exchange(arg, exec_flag);
    }
}

/*******************************************************************************
 * Routine to force download halos, if marked as dirty
*******************************************************************************/
void opp_mpi_force_halo_update_if_dirty(opp_set set, std::vector<opp_dat> dats, DeviceType device) {
    
    const int nargs = (int)dats.size();
    std::vector<opp_arg> args(nargs);

    for (int i = 0; i < nargs; i++) {
        args[i] = opp_arg_dat(dats[i], OPP_READ);
        args[i].idx = 2; // HACK to forcefully make halos to download
    }

    opp_mpi_halo_exchanges_grouped(set, nargs, args.data(), device);
    opp_mpi_halo_wait_all(nargs, args.data());
}

/*******************************************************************************
 * Routine to wait for all the MPI halo exchanges to complete
*******************************************************************************/
void opp_mpi_halo_wait_all(int nargs, opp_arg *args)
{
    __opp_mpi_host_halo_wait_all(nargs, args);
}