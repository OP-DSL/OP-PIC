
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

#include "opp_mpi.h"
#include "opp_hdf5.h"

#include "opp_hdf5_common.cpp"

// Use version 2 of H5Dopen H5Acreate and H5Dcreate
#define H5Dopen_vers 2
#define H5Acreate_vers 2
#define H5Dcreate_vers 2

MPI_Comm OPP_MPI_HDF5_WORLD;

int compute_local_size_weight(int global_size, int mpi_comm_size, int mpi_rank) {
    
    // TODO : refactor: Remove OPP_hybrid_gpu and OPP_hybrid_balance
    int OPP_hybrid_gpu = 1;
    int *hybrid_flags = (int *)malloc(mpi_comm_size * sizeof(int));
    MPI_Allgather(&OPP_hybrid_gpu, 1, MPI_INT, hybrid_flags, 1, MPI_INT,
                    OPP_MPI_HDF5_WORLD);
    double total = 0;
    for (int i = 0; i < mpi_comm_size; i++)
        total += hybrid_flags[i] == 1 ? OPP_hybrid_balance : 1.0;
    double cumulative = 0;
    for (int i = 0; i < mpi_rank; i++)
        cumulative += hybrid_flags[i] == 1 ? OPP_hybrid_balance : 1.0;

    int local_start = ((double)global_size) * (cumulative / total);
    int local_end = ((double)global_size) *
        ((cumulative + (OPP_hybrid_gpu ? OPP_hybrid_balance : 1.0)) / total);
    if (mpi_rank + 1 == mpi_comm_size)
        local_end = global_size; // make sure we don't have rounding problems
    int local_size = local_end - local_start;

    return local_size;
}

/*******************************************************************************
* Routine to read an opp_set from an hdf5 file
*******************************************************************************/
int opp_decl_set_hdf5_data(char const *file, char const *name) {
    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;  // file identifier
    hid_t plist_id; // property list identifier
    hid_t dset_id;  // dataset identifier

    if (file_exist(file) == 0) {
        printf("File %s does not exist .... aborting op_decl_set_hdf5()\n", file);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

    file_id = H5Fopen(file, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    // Create the dataset with default properties and close dataspace.
    dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        printf("Could not open dataset '%s' in file '%s'\n", name, file);
        H5Fclose(file_id);
        return -1;
    }

    // Create property list for collective dataset read.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    int g_size = 0;
    // read data
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, &g_size);

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);

    // calculate local size of set for this mpi process
    int l_size = compute_local_size_weight(g_size, comm_size, my_rank);
    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);

    return l_size;
}

opp_set opp_decl_set_hdf5(char const *file, char const *name) { 

    int l_size = 0;
    if (OPP_IS_VALID_PROCESS)
        l_size = opp_decl_set_hdf5_data(file, name);
    return opp_decl_set(l_size, name);
}

opp_set opp_decl_particle_set_hdf5(char const *file, char const *name, opp_set cells_set) { 

    int l_size = 0;
    if (OPP_IS_VALID_PROCESS)
        l_size = opp_decl_set_hdf5_data(file, name);
    return opp_decl_particle_set(l_size, name, cells_set);
}

// opp_set op_decl_set_hdf5_infer_size(char const *file, char const *name, char const *set_dataset_name) {
//   // printf("op_decl_set_hdf5_infer_size() called in op_mpi_hdf5.c\n");
//   // create new communicator
//   int my_rank, comm_size;
//   MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
//   MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
//   MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

//   // MPI variables
//   MPI_Info info = MPI_INFO_NULL;

//   // HDF5 APIs definitions
//   hid_t file_id;   // file identifier
//   hid_t plist_id;  // property list identifier
//   hid_t dset_id;   // dataset identifier
//   hid_t dataspace; // data space identifier
//   herr_t status;

//   if (file_exist(file) == 0) {
//     printf("File %s does not exist .... aborting op_decl_set_hdf5()\n",
//               file);
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }

//   // Set up file access property list with parallel I/O access
//   plist_id = H5Pcreate(H5P_FILE_ACCESS);
//   if (plist_id < 0) {
//     printf("H5Pcreate() failed\n");
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }
//   if (H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info) < 0) {
//     printf("H5Pset_fapl_mpio() failed\n");
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }

//   file_id = H5Fopen(file, H5F_ACC_RDONLY, plist_id);
//   H5Pclose(plist_id);
//   if (file_id < 0) {
//     printf("Could not obtain read access to file '%s'\n", file);
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }

//   if (!H5Lexists(file_id, set_dataset_name, H5P_DEFAULT)) {
//     printf("Dataset '%s' not found in file '%s'\n", set_dataset_name, file);
//     H5Fclose(file_id);
//     return NULL;
//   }

//   // Create the dataset with default properties and close dataspace.
//   // printf("Rank %d reading set %s\n", my_rank, name);
//   dset_id = H5Dopen(file_id, set_dataset_name, H5P_DEFAULT);
//   if (dset_id < 0) {
//     printf("Could not open dataset '%s' in file '%s'\n", set_dataset_name, file);
//     H5Fclose(file_id);
//     return NULL;
//   }

//   // Get size of dataset:
//   op_hdf5_dataset_properties dset_props;
//   status = get_dataset_properties(dset_id, &dset_props);
//   if (status < 0) {
//     printf("Could not access dataset '%s' in file '%s'\n", set_dataset_name, file);
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }
//   int g_size = dset_props.size;

//   H5Dclose(dset_id);
//   H5Fclose(file_id);

//   free((char*)dset_props.type_str);

//   // calculate local size of set for this mpi process
//   int l_size = compute_local_size_weight(g_size, comm_size, my_rank);
//   MPI_Comm_free(&OPP_MPI_HDF5_WORLD);

//   return op_decl_set(l_size, name);
// }

/*******************************************************************************
* Routine to read an opp_map from an hdf5 file
*******************************************************************************/

opp_map opp_decl_map_hdf5(opp_set from, opp_set to, int dim, char const *file,
                        char const *name) { OPP_RETURN_NULL_IF_INVALID_PROCESS;
  
    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;   // file identifier
    hid_t plist_id;  // property list identifier
    hid_t dset_id;   // dataset identifier
    hid_t dataspace; // data space identifier
    hid_t memspace;  // memory space identifier
    herr_t status;

    hsize_t count[2]; // hyperslab selection parameters
    hsize_t offset[2];

    if (file_exist(file) == 0) {
        printf("File %s does not exist .... aborting op_decl_map_hdf5()\n",
                file);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);
    file_id = H5Fopen(file, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5error_off(&old_func, &old_client_data);

    /*open data set*/
    dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        printf("opp_map with name : %s not found in file : %s \n", name, file);
        H5Fclose(file_id);
        return NULL;
    }

    op_hdf5_dataset_properties dset_props;
    status = get_dataset_properties(dset_id, &dset_props);
    if (status < 0) {
        printf("Could not get properties of dataset '%s' in file '%s'\n", name, file);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    int g_size = dset_props.size;

    // calculate local size of set for this mpi process
    int l_size = compute_local_size_weight(g_size, comm_size, my_rank);
    // check if size is accurate
    if (from->size != l_size) {
        printf(
            "map from set size %d in file %s and size %d do not match on rank %d\n",
            l_size, file, from->size, my_rank);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    int map_dim = dset_props.dim;
    if (map_dim != dim) {
        printf("map.dim %d in file %s and dim %d do not match\n", map_dim, file,
                dim);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    const char *typ = dset_props.type_str;

    // Restore previous error handler .. report hdf5 error stack automatically
    H5error_on(old_func, old_client_data);

    /*read in map in hyperslabs*/

    // Each process defines dataset in memory and reads from a hyperslab in the
    // file.
    int disp = 0;
    int *sizes = (int *)malloc(sizeof(int) * comm_size);
    MPI_Allgather(&l_size, 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
    for (int i = 0; i < my_rank; i++)
        disp = disp + sizes[i];
    free(sizes);

    count[0] = l_size;
    count[1] = dim;
    offset[0] = disp;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    // Select hyperslab in the file.
    dataspace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // initialize data buffer and read data
    int *map = 0;
    if (strcmp(typ, "int") == 0 || strcmp(typ, "integer(4)") == 0) {
        map = (int *)malloc(sizeof(int) * l_size * dim);
        H5Dread(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, map);
    } else if (strcmp(typ, "long") == 0) {
        map = (int *)malloc(sizeof(long) * l_size * dim);
        H5Dread(dset_id, H5T_NATIVE_LONG, memspace, dataspace, plist_id, map);
    } else if (strcmp(typ, "long long") == 0) {
        map = (int *)malloc(sizeof(long) * l_size * dim);
        H5Dread(dset_id, H5T_NATIVE_LLONG, memspace, dataspace, plist_id, map);
    } else {
        printf("unknown type\n");
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);

    free((char*)dset_props.type_str);

    opp_map new_map = opp_decl_map(from, to, dim, map, name);
    free(map);
    new_map->user_managed = 0;
    
    return new_map;
}

/*******************************************************************************
* Routine to read an opp_dat from an hdf5 file
*******************************************************************************/

char* opp_get_hdf5_dat_data(opp_set set, int dim, char const *type, char const *file, char const *name) {
  
    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;   // file identifier
    hid_t plist_id;  // property list identifier
    hid_t dset_id;   // dataset identifier
    hid_t dataspace; // data space identifier
    hid_t memspace;  // memory space identifier

    hsize_t count[2]; // hyperslab selection parameters
    hsize_t offset[2];
    // hid_t attr; // attribute identifier
    herr_t status;

    if (file_exist(file) == 0) {
        printf("File %s does not exist .... aborting op_decl_dat_hdf5()\n",
                file);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

    file_id = H5Fopen(file, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5error_off(&old_func, &old_client_data);

    /*open data set*/
    dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        if (!set->is_particle)
            printf("opp_dat with name : %s not found in file : %s \n", name, file);
        H5Fclose(file_id);
        return NULL;
    }

    op_hdf5_dataset_properties dset_props;
    status = get_dataset_properties(dset_id, &dset_props);
    if (status < 0) {
        printf("Could not get properties of dataset '%s' in file '%s'\n", name, file);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }
    //size_t type_size;

    int dat_dim = dset_props.dim;
    if (dat_dim != dim) {
        printf("dat.dim %d in file %s and dim %d do not match\n", dat_dim, file, dim);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    const char *typ = dset_props.type_str;
    if (!opp_type_equivalence(typ, type)) {
        if (OPP_DBG) 
            printf("dat.type %s in file %s and type %s do not match, performing automatic type conversion\n", 
                typ, file, type);
    }

    // Restore previous error handler .. report hdf5 error stack automatically
    H5error_on(old_func, old_client_data);

    /*read in dat in hyperslabs*/

    // Create the dataset with default properties and close dataspace.
    // Each process defines dataset in memory and reads from a hyperslab in the
    // file.
    int disp = 0;
    int *sizes = (int *)malloc(sizeof(int) * comm_size);
    MPI_Allgather(&(set->size), 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
    for (int i = 0; i < my_rank; i++)
        disp = disp + sizes[i];
    free(sizes);

    count[0] = set->size;
    count[1] = dim;
    offset[0] = disp;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    // Select hyperslab in the file.
    dataspace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // initialize data buffer and read in data
    char *data = 0;
    if (strcmp(type, "double") == 0 || strcmp(type, "double:soa") == 0 ||
        strcmp(type, "double precision") == 0 || strcmp(type, "real(8)") == 0) {
        data = (char *)malloc(set->size * dim * sizeof(double));
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, data);
        //type_size = sizeof(double);
    } else if (strcmp(type, "float") == 0 || strcmp(type, "float:soa") == 0 ||
                strcmp(type, "real(4)") == 0 || strcmp(type, "real") == 0) {
        data = (char *)malloc(set->size * dim * sizeof(float));
        H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id, data);
        //type_size = sizeof(float);

    } else if (strcmp(type, "int") == 0 || strcmp(type, "int:soa") == 0 ||
                strcmp(type, "int(4)") == 0 || strcmp(type, "integer") == 0 ||
                strcmp(type, "integer(4)") == 0) {
        data = (char *)malloc(set->size * dim * sizeof(int));
        H5Dread(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, data);
        //type_size = sizeof(int);
    } else {
        printf("unknown type\n");
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);

    free((char*)dset_props.type_str);
    
    return data;
}

opp_dat opp_decl_dat_hdf5(opp_set set, int dim, opp_data_type dtype, char const *file, char const *name) { OPP_RETURN_NULL_IF_INVALID_PROCESS;
    
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* data = opp_get_hdf5_dat_data(set, dim, type.c_str(), file, name);

    opp_dat new_dat = nullptr;
    new_dat = opp_decl_dat(set, dim, dtype, data, name);

    free(data);
    
    return new_dat;
}

/*******************************************************************************
* Routine to read in a constant from a named hdf5 file
*******************************************************************************/
void opp_get_const_hdf5(char const *name, int dim, char const *type, char *const_data, 
                            char const *file_name) { OPP_RETURN_IF_INVALID_PROCESS;
    
    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;   // file identifier
    hid_t plist_id;  // property list identifier
    hid_t dset_id;   // dataset identifier
    // hid_t dataspace; // data space identifier
    // hid_t attr;      // attribute identifier
    hid_t status;

    if (file_exist(file_name) == 0) {
        printf("File %s does not exist .... aborting op_get_const_hdf5()\n", file_name);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

    file_id = H5Fopen(file_name, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    // find dimension of this constant with available attributes
    int const_dim = 0;
    dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    if (dset_id < 0) {
        printf("dataset with '%s' not found in file '%s' \n", name, file_name);
        H5Fclose(file_id);
        const_data = NULL;
        return;
    }

    // Create property list for collective dataset read.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    op_hdf5_dataset_properties dset_props;
    status = get_dataset_properties(dset_id, &dset_props);
    if (status < 0) {
        printf("Could not get properties of dataset '%s' in file '%s'\n", name, file_name);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    const_dim = dset_props.size;
    if (const_dim != dim) {
        printf("dim of constant %d in file %s and requested dim %d do not match\n",
                const_dim, file_name, dim);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    const char* typ = dset_props.type_str;
    if (!opp_type_equivalence(typ, type)) {
        if (OPP_DBG) printf(
            "type of constant %s in file %s and requested type %s do not match, performing automatic type conversion\n",
            typ, file_name, type);
        typ = type;
    }
    H5Dclose(dset_id);

    // Create the dataset with default properties and close dataspace.
    dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
    // dataspace = H5Dget_space(dset_id);

    char *data = nullptr;
    // initialize data buffer and read data
    if (strcmp(typ, "int") == 0 || strcmp(typ, "int(4)") == 0 ||
        strcmp(typ, "integer") == 0 || strcmp(typ, "integer(4)") == 0) {
        data = (char *)malloc(sizeof(int) * const_dim);
        H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, data);
        memcpy((void *)const_data, (void *)data, sizeof(int) * const_dim);
    } else if (strcmp(typ, "long") == 0) {
        data = (char *)malloc(sizeof(long) * const_dim);
        H5Dread(dset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, plist_id, data);
        memcpy((void *)const_data, (void *)data, sizeof(long) * const_dim);
    } else if (strcmp(typ, "long long") == 0) {
        data = (char *)malloc(sizeof(long long) * const_dim);
        H5Dread(dset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, plist_id, data);
        memcpy((void *)const_data, (void *)data, sizeof(long long) * const_dim);
    } else if (strcmp(typ, "float") == 0 || strcmp(typ, "real(4)") == 0 ||
                strcmp(typ, "real") == 0) {
        data = (char *)malloc(sizeof(float) * const_dim);
        H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, plist_id, data);
        memcpy((void *)const_data, (void *)data, sizeof(float) * const_dim);
    } else if (strcmp(typ, "double") == 0 ||
                strcmp(typ, "double precision") == 0 ||
                strcmp(typ, "real(8)") == 0) {
        data = (char *)malloc(sizeof(double) * const_dim);
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, data);
        memcpy((void *)const_data, (void *)data, sizeof(double) * const_dim);
    } else {
        printf("Unknown type in file %s for constant %s\n", file_name, name);
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    free(data);

    free((char*)dset_props.type_str);

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);
}

/*******************************************************************************
* Routine to write all to a named hdf5 file
*******************************************************************************/
void opp_dump_to_hdf5(char const *file_name) { OPP_RETURN_IF_INVALID_PROCESS;

    opp_printf("opp_dump_to_hdf5", "Writing to %s", file_name);

    opp_profiler->start("Dump_hdf5");

    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;     // file identifier
    hid_t plist_id;    // property list identifier
    hid_t dset_id = 0; // dataset identifier
    hid_t dataspace;   // data space identifier
    hid_t memspace;    // memory space identifier

    hsize_t dimsf[2]; // dataset dimensions
    hsize_t count[2]; // hyperslab selection parameters
    hsize_t offset[2];

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

    // Create a new file collectively and release property list identifier.
    file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*loop over all the opp_sets and write them to file*/
    for (int s = 0; s < (int)opp_sets.size(); s++) {
        
        opp_set set = opp_sets[s];

        // Create the dataspace for the dataset.
        hsize_t dimsf_set[] = {1};
        dataspace = H5Screate_simple(1, dimsf_set, NULL);

        // Create the dataset with default properties and close dataspace.
        dset_id = H5Dcreate(file_id, set->name, H5T_NATIVE_INT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Create property list for collective dataset write.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        int size = 0;
        int *sizes = (int *)malloc(sizeof(int) * comm_size);
        MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
        for (int i = 0; i < comm_size; i++)
            size = size + sizes[i];

        // write data
        H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, &size);
        H5Sclose(dataspace);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    /*loop over all the opp_maps and write them to file*/
    for (int m = 0; m < (int)opp_maps.size(); m++) {
        
        opp_map map = opp_maps[m];

        if (map->dim == 0)
            continue;

        // find total size of map
        int *sizes = (int *)malloc(sizeof(int) * comm_size);
        int g_size = 0;
        MPI_Allgather(&map->from->size, 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
        for (int i = 0; i < comm_size; i++)
            g_size = g_size + sizes[i];

        if (g_size == 0)
            continue;

        // Create the dataspace for the dataset.
        dimsf[0] = g_size;
        dimsf[1] = map->dim;
        dataspace = H5Screate_simple(2, dimsf, NULL);

        // Each process defines dataset in memory and writes it to a hyperslab in the file.
        int disp = 0;
        for (int i = 0; i < my_rank; i++)
            disp = disp + sizes[i];
        count[0] = map->from->size;
        count[1] = dimsf[1];
        offset[0] = disp;
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        // Select hyperslab in the file.
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Create property list for collective dataset write.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // create map path
        create_path(map->name, file_id);

        // Create the dataset with default properties and close dataspace.
        if (sizeof(map->map[0]) == sizeof(int)) {
                
            dset_id = H5Dcreate(file_id, map->name, H5T_NATIVE_INT, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, map->map);

        } else if (sizeof(map->map[0]) == sizeof(long)) {

            dset_id = H5Dcreate(file_id, map->name, H5T_NATIVE_LONG, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_LONG, memspace, dataspace, plist_id, map->map);

        } else if (sizeof(map->map[0]) == sizeof(long long)) {
            dset_id = H5Dcreate(file_id, map->name, H5T_NATIVE_LLONG, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_LLONG, memspace, dataspace, plist_id, map->map);
        }

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        free(sizes);

        /*attach attributes to map*/

        // open existing data set
        dset_id = H5Dopen(file_id, map->name, H5P_DEFAULT);
        // create the data space for the attribute
        hsize_t dims = 1;
        dataspace = H5Screate_simple(1, &dims, NULL);

        // Create an int attribute - size
        hid_t attribute = H5Acreate(dset_id, "size", H5T_NATIVE_INT, dataspace,
                                    H5P_DEFAULT, H5P_DEFAULT);
        // Write the attribute data.
        H5Awrite(attribute, H5T_NATIVE_INT, &g_size);
        // Close the attribute.
        H5Aclose(attribute);

        // Create an int attribute - dimension
        attribute = H5Acreate(dset_id, "dim", H5T_NATIVE_INT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT);
        // Write the attribute data.
        H5Awrite(attribute, H5T_NATIVE_INT, &map->dim);
        // Close the attribute.
        H5Aclose(attribute);
        H5Sclose(dataspace);

        // Create a string attribute - type
        dataspace = H5Screate(H5S_SCALAR);
        hid_t atype = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, 10);
        attribute = H5Acreate(dset_id, "type", atype, dataspace, H5P_DEFAULT, H5P_DEFAULT);

        if (sizeof(map->map[0]) == sizeof(int))
            H5Awrite(attribute, atype, "int");
        if (sizeof(map->map[0]) == sizeof(long))
            H5Awrite(attribute, atype, "long");
        if (sizeof(map->map[0]) == sizeof(long long))
            H5Awrite(attribute, atype, "long long");

        H5Aclose(attribute);
        // Close the dataspace
        H5Sclose(dataspace);
        // Close to the dataset.
        H5Dclose(dset_id);
    }

    /*loop over all the opp_dats and write them to file*/
    for (int k = 0; k < (int)opp_dats.size(); k++) {
        
        opp_dat dat = opp_dats[k];
        
        if (dat->size == 0)
            continue;
        
        // avoid temp dats
        if (strncmp(dat->name, "AUTO_DAT_", 8) == 0)
            continue;

        // find total size of dat
        int *sizes = (int *)malloc(sizeof(int) * comm_size);
        int g_size = 0;
        MPI_Allgather(&dat->set->size, 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
        for (int i = 0; i < comm_size; i++)
            g_size = g_size + sizes[i];
        
        if (g_size == 0)
            continue;

        // Create the dataspace for the dataset.
        dimsf[0] = g_size;
        dimsf[1] = dat->dim;
        dataspace = H5Screate_simple(2, dimsf, NULL);

        // Each process defines dataset in memory and writes it to a hyperslab in the file.
        int disp = 0;
        for (int i = 0; i < my_rank; i++)
            disp = disp + sizes[i];
        count[0] = dat->set->size;
        count[1] = dimsf[1];
        offset[0] = disp;
        offset[1] = 0;
        memspace = H5Screate_simple(2, count, NULL);

        // Select hyperslab in the file.
        H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Create property list for collective dataset write.
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // create dateset path
        create_path(dat->name, file_id);

        // Create the dataset with default properties and close dataspace.
        if (strcmp(dat->type, "double") == 0 ||
            strcmp(dat->type, "double:soa") == 0 ||
            strcmp(dat->type, "double precision") == 0 ||
            strcmp(dat->type, "real(8)") == 0) {
            
            dset_id = H5Dcreate(file_id, dat->name, H5T_NATIVE_DOUBLE, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, dat->data);

        } else if (strcmp(dat->type, "float") == 0 ||
                strcmp(dat->type, "float:soa") == 0 ||
                strcmp(dat->type, "real(4)") == 0 ||
                strcmp(dat->type, "real") == 0) {

            dset_id = H5Dcreate(file_id, dat->name, H5T_NATIVE_FLOAT, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id, dat->data);

        } else if (strcmp(dat->type, "int") == 0 ||
                strcmp(dat->type, "int:soa") == 0 ||
                strcmp(dat->type, "int(4)") == 0 ||
                strcmp(dat->type, "integer") == 0 ||
                strcmp(dat->type, "integer(4)") == 0) {

            dset_id = H5Dcreate(file_id, dat->name, H5T_NATIVE_INT, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, dat->data);

        } else if ((strcmp(dat->type, "long") == 0) ||
                (strcmp(dat->type, "long:soa") == 0)) {

            dset_id = H5Dcreate(file_id, dat->name, H5T_NATIVE_LONG, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_LONG, memspace, dataspace, plist_id, dat->data);

        } else if ((strcmp(dat->type, "long long") == 0) ||
                (strcmp(dat->type, "long long:soa") == 0)) {

            dset_id = H5Dcreate(file_id, dat->name, H5T_NATIVE_LLONG, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_LLONG, memspace, dataspace, plist_id, dat->data);

        } else {

            printf("Unknown type - in op_dump_to_hdf5() writing opp_dats\n");
            MPI_Abort(OPP_MPI_HDF5_WORLD, 2);

        }

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(dataspace);
        free(sizes);

        /*attach attributes to dat*/

        // open existing data set
        dset_id = H5Dopen(file_id, dat->name, H5P_DEFAULT);
        
        // create the data space for the attribute
        hsize_t dims = 1;
        dataspace = H5Screate_simple(1, &dims, NULL);

        // Create an int attribute - size
        hid_t attribute = H5Acreate(dset_id, "size", H5T_NATIVE_INT, dataspace,
                                    H5P_DEFAULT, H5P_DEFAULT);
        // Write the attribute data.
        H5Awrite(attribute, H5T_NATIVE_INT, &dat->size);
        
        // Close the attribute.
        H5Aclose(attribute);

        // Create an int attribute - dimension
        attribute = H5Acreate(dset_id, "dim", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        
        // Write the attribute data.
        H5Awrite(attribute, H5T_NATIVE_INT, &dat->dim);
        H5Aclose(attribute);
        H5Sclose(dataspace);

        // Create a string attribute - type
        dataspace = H5Screate(H5S_SCALAR);
        hid_t atype = H5Tcopy(H5T_C_S1);
        int attlen = strlen(dat->type);
        H5Tset_size(atype, attlen);

        attribute = H5Acreate(dset_id, "type", atype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute, atype, dat->type);

        H5Aclose(attribute);
        H5Sclose(dataspace);
        H5Dclose(dset_id);
    }
    H5Fclose(file_id);

    opp_profiler->end("Dump_hdf5");

    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);
}

// /*******************************************************************************
// * Routine to write a constant to a named hdf5 file
// *******************************************************************************/
// void op_write_const_hdf5(char const *name, int dim, char const *type,
//                          char *const_data, char const *file_name) {
//   // letting know that writing is happening ...
//   printf("Writing '%s' to file '%s'\n", name, file_name);

//   // create new communicator
//   int my_rank, comm_size;
//   MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
//   MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
//   MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

//   // MPI variables
//   MPI_Info info = MPI_INFO_NULL;

//   // HDF5 APIs definitions
//   hid_t file_id;   // file identifier
//   hid_t dset_id;   // dataset identifier
//   hid_t plist_id;  // property list identifier
//   hid_t dataspace; // data space identifier

//   // Set up file access property list with parallel I/O access
//   plist_id = H5Pcreate(H5P_FILE_ACCESS);
//   H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

//   if (file_exist(file_name) == 0) {
//     if (OPP_DBG) {
//       printf("File %s does not exist .... creating file\n", file_name);
//     }
//     file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
//     H5Fclose(file_id);
//   }

//   /* Open the existing file. */
//   file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
//   H5Pclose(plist_id);

//   // Create the dataspace for the dataset.
//   hsize_t dims_of_const = {(hsize_t)dim};
//   dataspace = H5Screate_simple(1, &dims_of_const, NULL);

//   // Create property list for collective dataset write.
//   plist_id = H5Pcreate(H5P_DATASET_XFER);
//   H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

//   // Create the dataset with default properties
//   if (strcmp(type, "double") == 0 || strcmp(type, "double:soa") == 0 ||
//       strcmp(type, "double precision") == 0 || strcmp(type, "real(8)") == 0) {
//     dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, dataspace,
//                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     // write data
//     H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace, plist_id,
//              const_data);
//     H5Dclose(dset_id);
//   } else if (strcmp(type, "float") == 0 || strcmp(type, "float:soa") == 0 ||
//              strcmp(type, "real(4)") == 0 || strcmp(type, "real") == 0) {
//     dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT,
//                         H5P_DEFAULT, H5P_DEFAULT);
//     // write data
//     H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, dataspace, plist_id,
//              const_data);
//     H5Dclose(dset_id);
//   } else if (strcmp(type, "int") == 0 || strcmp(type, "int:soa") == 0 ||
//              strcmp(type, "int(4)") == 0 || strcmp(type, "integer") == 0 ||
//              strcmp(type, "integer(4)") == 0) {
//     dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
//                         H5P_DEFAULT, H5P_DEFAULT);
//     // write data
//     H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, dataspace, plist_id, const_data);
//     H5Dclose(dset_id);
//   } else if ((strcmp(type, "long") == 0) || (strcmp(type, "long:soa") == 0)) {
//     dset_id = H5Dcreate(file_id, name, H5T_NATIVE_LONG, dataspace, H5P_DEFAULT,
//                         H5P_DEFAULT, H5P_DEFAULT);
//     // write data
//     H5Dwrite(dset_id, H5T_NATIVE_LONG, H5S_ALL, dataspace, plist_id,
//              const_data);
//     H5Dclose(dset_id);
//   } else if ((strcmp(type, "long long") == 0) ||
//              (strcmp(type, "long long:soa") == 0)) {
//     dset_id = H5Dcreate(file_id, name, H5T_NATIVE_LLONG, dataspace, H5P_DEFAULT,
//                         H5P_DEFAULT, H5P_DEFAULT);
//     // write data
//     H5Dwrite(dset_id, H5T_NATIVE_LLONG, H5S_ALL, dataspace, plist_id,
//              const_data);
//     H5Dclose(dset_id);
//   } else {
//     printf("Unknown type %s for constant %s: cannot write constant to file\n",
//               type, name);
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }

//   H5Pclose(plist_id);
//   H5Sclose(dataspace);

//   /*attach attributes to constant*/

//   // open existing data set
//   dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
//   // create the data space for the attribute
//   dims_of_const = 1;
//   dataspace = H5Screate_simple(1, &dims_of_const, NULL);

//   // Create an int attribute - dimension
//   hid_t attribute = H5Acreate(dset_id, "dim", H5T_NATIVE_INT, dataspace,
//                               H5P_DEFAULT, H5P_DEFAULT);
//   // Write the attribute data.
//   H5Awrite(attribute, H5T_NATIVE_INT, &dim);
//   // Close the attribute.
//   H5Aclose(attribute);
//   H5Sclose(dataspace);

//   // Create a string attribute - type
//   dataspace = H5Screate(H5S_SCALAR);
//   hid_t atype = H5Tcopy(H5T_C_S1);

//   int attlen = strlen(type);
//   H5Tset_size(atype, attlen);
//   attribute =
//       H5Acreate(dset_id, "type", atype, dataspace, H5P_DEFAULT, H5P_DEFAULT);

//   if (strcmp(type, "double") == 0 || strcmp(type, "double:soa") == 0 ||
//       strcmp(type, "double precision") == 0 || strcmp(type, "real(8)") == 0)
//     H5Awrite(attribute, atype, "double");
//   else if (strcmp(type, "int") == 0 || strcmp(type, "int:soa") == 0 ||
//            strcmp(type, "int(4)") == 0 || strcmp(type, "integer") == 0 ||
//            strcmp(type, "integer(4)") == 0)
//     H5Awrite(attribute, atype, "int");
//   else if (strcmp(type, "long") == 0)
//     H5Awrite(attribute, atype, "long");
//   else if (strcmp(type, "long long") == 0)
//     H5Awrite(attribute, atype, "long long");
//   else if (strcmp(type, "float") == 0 || strcmp(type, "float:soa") == 0 ||
//            strcmp(type, "real(4)") == 0 || strcmp(type, "real") == 0)
//     H5Awrite(attribute, atype, "float");
//   else {
//     printf("Unknown type %s for constant %s: cannot write constant to file\n",
//               type, name);
//     MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
//   }

//   H5Aclose(attribute);
//   H5Sclose(dataspace);
//   H5Dclose(dset_id);

//   H5Fclose(file_id);
//   MPI_Comm_free(&OPP_MPI_HDF5_WORLD);
// }

/*******************************************************************************
* Routine to write an opp_dat to a named hdf5 file,
* if file does not exist, creates it
* if the data set given at path does not exists in file creates data set
*******************************************************************************/

void opp_fetch_data_hdf5(opp_dat data, char const *file_name, char const *path_name) { OPP_RETURN_IF_INVALID_PROCESS;

    // letting know that writing is happening ...
    opp_printf("opp_fetch_data_hdf5", "Fetching from '%s' to '%s'\n", path_name, file_name);

    // fetch data based on the backend
    opp_dat dat = opp_fetch_data(data);

    // create new communicator
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_HDF5_WORLD);
    MPI_Comm_rank(OPP_MPI_HDF5_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_HDF5_WORLD, &comm_size);

    // MPI variables
    MPI_Info info = MPI_INFO_NULL;

    // HDF5 APIs definitions
    hid_t file_id;     // file identifier
    hid_t dset_id = 0; // dataset identifier
    hid_t dataspace;   // data space identifier
    hid_t plist_id;    // property list identifier
    hid_t memspace;    // memory space identifier
    // hid_t attr;        // attribute identifier

    hsize_t dimsf[2]; // dataset dimensions
    hsize_t count[2]; // hyperslab selection parameters
    hsize_t offset[2];

    H5E_auto_t old_func;
    void *old_client_data;

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, OPP_MPI_HDF5_WORLD, info);

    if (file_exist(file_name) == 0) {

        MPI_Barrier(OPP_MPI_WORLD);
        if (OPP_DBG) {
            printf("File %s does not exist .... creating file\n", file_name);
        }
        MPI_Barrier(OPP_MPI_HDF5_WORLD);
        if (OPP_rank == OPP_ROOT) {
            FILE *fp;
            fp = fopen(file_name, "w");
            fclose(fp);
        }
        MPI_Barrier(OPP_MPI_HDF5_WORLD);
        file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    
    } 
    else {

        if (OPP_DBG) {
            printf("File %s exists .... checking for dataset %s in file\n",
                    file_name, path_name);
        }
        file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);

        H5error_off(&old_func, &old_client_data);
        herr_t status = H5Gget_objinfo(file_id, path_name, 0, NULL);
        H5error_on(old_func, old_client_data);

        if (status == 0) {
            if (OPP_DBG) {
                printf("opp_dat %s exists in the file ... updating data\n", path_name);
            }

            dset_id = H5Dopen(file_id, path_name, H5P_DEFAULT);

            op_hdf5_dataset_properties dset_props;
            herr_t status = get_dataset_properties(dset_id, &dset_props);
            if (status < 0) {
                printf("Could not get properties of dataset '%s' in file '%s'\n",
                        path_name, file_name);
                MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
            }

            // find element size of this dat with available attributes
            size_t dat_size = 0;
            dat_size = dset_props.elem_bytes;
            if (dat_size != (size_t)dat->size) {
                printf("dat.size %zu in file %s and dim %d do not match ... aborting\n",
                    dat_size, file_name, dat->dim);
                MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
            }

            // find dim with available attributes
            int dat_dim = 0;
            dat_dim = dset_props.dim;
            if (dat_dim != dat->dim) {
                printf("dat.dim %d in file %s and dim %d do not match ... aborting\n",
                        dat_dim, file_name, dat->dim);
                MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
            }

            // find type with available attributes
            dataspace = H5Screate(H5S_SCALAR);
            const char* typ = dset_props.type_str;
            char typ_soa[50];
            sprintf(typ_soa, "%s:soa", typ);
            if (!opp_type_equivalence(typ, dat->type)) {
                printf("dat.type %s in file %s and type %s do not match\n", typ,
                        file_name, dat->type);
                MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
            }

            //
            // all good .. we can update existing dat now
            //

            // find total size of dat
            int *sizes = (int *)malloc(sizeof(int) * comm_size);
            int g_size = 0;
            MPI_Allgather(&dat->set->size, 1, MPI_INT, sizes, 1, MPI_INT, OPP_MPI_HDF5_WORLD);
            for (int i = 0; i < comm_size; i++)
                g_size = g_size + sizes[i];

            // Create the dataspace for the dataset.
            dimsf[0] = g_size;
            dimsf[1] = dat->dim;

            // Each process defines dataset in memory and writes it to a hyperslab
            // in the file.
            int disp = 0;
            for (int i = 0; i < my_rank; i++)
                disp = disp + sizes[i];
            count[0] = dat->set->size;
            count[1] = dimsf[1];
            offset[0] = disp;
            offset[1] = 0;
            memspace = H5Screate_simple(2, count, NULL);

            // Select hyperslab in the file.
            dataspace = H5Dget_space(dset_id);
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

            // Create property list for collective dataset write.
            H5Pclose(plist_id);
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

            // write data
            if (strcmp(dat->type, "double") == 0 ||
                strcmp(dat->type, "double:soa") == 0 ||
                strcmp(dat->type, "double precision") == 0 ||
                strcmp(dat->type, "real(8)") == 0)

                H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, dat->data);

            else if (strcmp(dat->type, "float") == 0 ||
                    strcmp(dat->type, "float:soa") == 0 ||
                    strcmp(dat->type, "real(4)") == 0 ||
                    strcmp(dat->type, "real") == 0)

                H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id, dat->data);

            else if (strcmp(dat->type, "int") == 0 ||
                    strcmp(dat->type, "int:soa") == 0 ||
                    strcmp(dat->type, "int(4)") == 0 ||
                    strcmp(dat->type, "integer") == 0 ||
                    strcmp(dat->type, "integer(4)") == 0)
                    
                H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, dat->data);

            else if ((strcmp(dat->type, "long") == 0) ||
                    (strcmp(dat->type, "long:soa") == 0))

                H5Dwrite(dset_id, H5T_NATIVE_LONG, memspace, dataspace, plist_id, dat->data);

            else if ((strcmp(dat->type, "long long") == 0) ||
                    (strcmp(dat->type, "long long:soa") == 0))

                H5Dwrite(dset_id, H5T_NATIVE_LLONG, memspace, dataspace, plist_id, dat->data);

            else {

                printf("Unknown type in op_fetch_data_hdf5_file()\n");
                MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
            }

            H5Dclose(dset_id);
            H5Pclose(plist_id);
            H5Sclose(memspace);
            H5Sclose(dataspace);
            H5Fclose(file_id);
            free(sizes);

            // free the temp opp_dat used for this write
            free(dat->data);
            free(dat->set);
            free(dat);

            MPI_Comm_free(&OPP_MPI_HDF5_WORLD);
            return;
        } 
        else {
            if (OPP_DBG) {
                printf("opp_dat %s does not exists in the file ... creating data set\n",
                        path_name);
            }
        }
    }

    //
    // new data set ...
    //

    H5Pclose(plist_id);

    // find total size of dat
    int *sizes = (int *)malloc(sizeof(int) * comm_size);
    int g_size = 0;
    MPI_Allgather(&dat->set->size, 1, MPI_INT, sizes, 1, MPI_INT,
                    OPP_MPI_HDF5_WORLD);
    for (int i = 0; i < comm_size; i++)
        g_size = g_size + sizes[i];

    // Create the dataspace for the dataset.
    dimsf[0] = g_size;
    dimsf[1] = dat->dim;
    dataspace = H5Screate_simple(2, dimsf, NULL);

    // Each process defines dataset in memory and writes it to a hyperslab in the
    // file.
    int disp = 0;
    for (int i = 0; i < my_rank; i++)
        disp = disp + sizes[i];
    count[0] = dat->set->size;
    count[1] = dimsf[1];
    offset[0] = disp;
    offset[1] = 0;
    memspace = H5Screate_simple(2, count, NULL);

    // Select hyperslab in the file.
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // create dataset path
    create_path(path_name, file_id);

    // Create the dataset with default properties and close dataspace.
    if (strcmp(dat->type, "double") == 0 ||
        strcmp(dat->type, "double:soa") == 0 ||
        strcmp(dat->type, "double precision") == 0 ||
        strcmp(dat->type, "real(8)") == 0) {
        dset_id = H5Dcreate(file_id, path_name, H5T_NATIVE_DOUBLE, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, dat->data);
    } else if (strcmp(dat->type, "float") == 0 ||
                strcmp(dat->type, "float:soa") == 0 ||
                strcmp(dat->type, "real(4)") == 0 ||
                strcmp(dat->type, "real") == 0) {
        dset_id = H5Dcreate(file_id, path_name, H5T_NATIVE_FLOAT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, dataspace, plist_id,
                dat->data);
    } else if (strcmp(dat->type, "int") == 0 ||
                strcmp(dat->type, "int:soa") == 0 ||
                strcmp(dat->type, "int(4)") == 0 ||
                strcmp(dat->type, "integer") == 0 ||
                strcmp(dat->type, "integer(4)") == 0) {
        dset_id = H5Dcreate(file_id, path_name, H5T_NATIVE_INT, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, dataspace, plist_id, dat->data);
    } else if ((strcmp(dat->type, "long") == 0) ||
                (strcmp(dat->type, "long:soa") == 0)) {
        dset_id = H5Dcreate(file_id, path_name, H5T_NATIVE_LONG, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_LONG, memspace, dataspace, plist_id,
                dat->data);
    } else if ((strcmp(dat->type, "long long") == 0) ||
                (strcmp(dat->type, "long long:soa") == 0)) {
        dset_id = H5Dcreate(file_id, path_name, H5T_NATIVE_LLONG, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_LLONG, memspace, dataspace, plist_id,
                dat->data);
    } else {
        printf("Unknown type in op_fetch_data_hdf5_file()\n");
        MPI_Abort(OPP_MPI_HDF5_WORLD, 2);
    }

    H5Dclose(dset_id);
    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(dataspace);
    free(sizes);

    /*attach attributes to dat*/

    // open existing data set
    dset_id = H5Dopen(file_id, path_name, H5P_DEFAULT);

    // create the data space for the attribute
    hsize_t dims = 1;
    dataspace = H5Screate_simple(1, &dims, NULL);

    // Create an int attribute - size
    hid_t attribute = H5Acreate(dset_id, "size", H5T_NATIVE_INT, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute, H5T_NATIVE_INT, &dat->size);
    H5Aclose(attribute);

    // Create an int attribute - dimension
    attribute = H5Acreate(dset_id, "dim", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute, H5T_NATIVE_INT, &dat->dim);
    H5Aclose(attribute);
    H5Sclose(dataspace);

    // Create an string attribute - type
    dataspace = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    int attlen = strlen(dat->type);
    H5Tset_size(atype, attlen);

    attribute = H5Acreate(dset_id, "type", atype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute, atype, dat->type);
    H5Aclose(attribute);

    H5Sclose(dataspace);
    H5Dclose(dset_id);
    H5Fclose(file_id);

    // free the temp opp_dat used for this write
    free(dat->data);
    free(dat->set);
    free(dat);

    MPI_Comm_free(&OPP_MPI_HDF5_WORLD);
}

/*******************************************************************************
* Routine to write an opp_dat to a named hdf5 file,
* if file does not exist, creates it
* if the data set does not exists in file creates data set
*******************************************************************************/

void opp_fetch_data_hdf5_file(opp_dat dat, char const *file_name) { OPP_RETURN_IF_INVALID_PROCESS;

    opp_fetch_data_hdf5(dat, file_name, dat->name);
}

// void op_fetch_data_hdf5_file_ptr(char *data, const char *file_name) {
//   opp_dat_entry *item;
//   opp_dat_entry *tmp_item;
//   opp_dat item_dat = NULL;
//   for (item = TAILQ_FIRST(&opp_dat_list); item != NULL; item = tmp_item) {
//     tmp_item = TAILQ_NEXT(item, entries);
//     // printf("Available opp_dat %s with pointer %p\n", item->dat->name,
//     // item->dat->data);
//     if (item->orig_ptr == data) {
//       // printf("%s(%p), ", item->dat->name, item->dat->data);
//       item_dat = item->dat;
//       break;
//     }
//   }
//   // printf("\n");
//   if (item_dat == NULL) {
//     printf("ERROR in op_partition: opp_dat not found for dat with %p pointer\n",
//            data);
//   }

//   op_fetch_data_hdf5(item_dat, file_name, item_dat->name);
// }

/*******************************************************************************
* Routine to write an opp_dat to a named hdf5 file, where dataset_name is a
* freely chosen path inside the h5-file of the dataset.
* If file does not exist, creates it
* If the data set does not exists in file creates data set
*******************************************************************************/

void opp_fetch_data_hdf5_file_path(opp_dat dat, char const *file_name, char const *path_name) { OPP_RETURN_IF_INVALID_PROCESS;

    opp_fetch_data_hdf5(dat, file_name, path_name);
}
