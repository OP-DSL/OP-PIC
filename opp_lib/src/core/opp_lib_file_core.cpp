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

#include <opp_lib_core.h>

//****************************************
void opp_print_dat_to_txtfile_core(opp_dat dat, const char *file_prefix, const char *file_suffix)
{ 
    opp_profiler->start("PrintFile");

    const std::string file_name = std::string("files/") + file_prefix + "_" + file_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) {
        opp_printf("opp_print_dat_to_txtfile_core", "can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d -- %d %d\n", 
            dat->set->size, dat->dim, dat->set->exec_size, dat->set->nonexec_size) < 0) {
        opp_printf("opp_print_dat_to_txtfile_core", "error writing to %s\n", file_name.c_str());
        exit(2);
    }

    const size_t row_count = (size_t)dat->set->size + dat->set->exec_size + dat->set->nonexec_size;
    for (size_t i = 0; i < row_count; i++) {
        for (int j = 0; j < dat->dim; j++) {
            if (strcmp(dat->type, "double") == 0) {
                if (((double *)dat->data)[i * dat->dim + j] == -0.0) { 
                    ((double *)dat->data)[i * dat->dim + j] = +0.0; 
                }

                if (fprintf(fp, " %2.25lE", ((double *)dat->data)[i * dat->dim + j]) < 0) {
                    printf("\topp_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if (strcmp(dat->type, "float") == 0) {
                if (fprintf(fp, " %f", ((float *)dat->data)[i * dat->dim + j]) < 0) {
                    printf("\topp_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if (strcmp(dat->type, "int") == 0) {
                if (fprintf(fp, " %d", ((int *)dat->data)[i * dat->dim + j]) < 0) {
                    printf("\topp_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else if (strcmp(dat->type, "long") == 0) {
                if (fprintf(fp, " %ld", ((long *)dat->data)[i * dat->dim + j]) < 0) {
                    printf("\topp_print_dat_to_txtfile_core error writing to %s\n", file_name.c_str());
                    exit(2);
                }
            } 
            else {
                printf("\topp_print_dat_to_txtfile_core Unknown type %s, cannot be written to file %s\n", 
                    dat->type, file_name.c_str());
                exit(2);
            }
        }

        fprintf(fp, "\n");

        if ((i + 1) == (size_t)dat->set->size) 
            fprintf(fp, "import_exec_below ****************************************\n");
        if ((i + 1) == (size_t)(dat->set->size + dat->set->exec_size)) 
            fprintf(fp, "import_non_exec_below ****************************************\n");
    }
    
    fclose(fp);

    opp_profiler->end("PrintFile");
}

//****************************************
void opp_dump_dat_core(opp_dat data) 
{
    fflush(stdout);

    if (data != NULL) {
        for (int i = 0; i < data->set->size; i++) {
            for (int j = 0; j < data->dim; j++) {
                if (strncmp("double", data->type, 6) == 0) {
                    printf("\t %+2.25lE", ((double *)data->data)[i * data->dim + j]);
                } 
                else if (strncmp("real", data->type, 4) == 0) {
                    printf("\t %f", ((float *)data->data)[i * data->dim + j]);
                } 
                else if (strncmp("integer", data->type, 7) == 0) {
                    printf("\t %d", data->data[i * data->dim + j]);
                } 
                else {
                    printf("\topp_dump_dat_core Unsupported type for dumping %s\n", data->type);
                    exit(0);
                }
            }

            printf("\t\n");
        }
    }

    fflush(stdout);
}

//****************************************
void opp_print_map_to_txtfile_core(opp_map map, const char *file_prefix, const char *file_suffix)
{ 
    if (map->from->is_particle) {
        opp_print_dat_to_txtfile_core(map->p2c_dat, file_prefix, file_suffix);
        return;
    }

    opp_profiler->start("PrintFile");

    const std::string file_name = std::string("files/") + file_prefix + "_" + file_suffix; 

    FILE *fp;
    if ((fp = fopen(file_name.c_str(), "w")) == NULL) {
        printf("\topp_print_map_to_txtfile_core can't open file %s\n", file_name.c_str());
        exit(2);
    }

    if (fprintf(fp, "%d %d -- %d %d\n", 
                map->from->size, map->dim, map->from->exec_size, map->from->nonexec_size) < 0) {
        printf("\topp_print_map_to_txtfile_core error writing to %s\n", file_name.c_str());
        exit(2);
    }

    for (int i = 0; i < map->from->size + map->from->exec_size; i++) {
        for (int j = 0; j < map->dim; j++) {
            if (fprintf(fp, " %d", ((int *)map->map)[i * map->dim + j]) < 0) {
                printf("\topp_print_map_to_txtfile_core error writing to %s\n", file_name.c_str());
                exit(2);
            }
        }

        fprintf(fp, "\n");

        if (i+1 == map->from->size) 
            fprintf(fp, "import_exec_below ****************************************\n");
        if (i+1 == map->from->size + map->from->exec_size) 
            fprintf(fp, "import_non_exec_below ****************************************\n");
    }
    
    fclose(fp);

    opp_profiler->end("PrintFile");
}
 
//****************************************
void* opp_load_from_file_core(const char* file_name, int set_size, int dim, char const *type, int size)
{
    int fsize = -1, fdim = -1;
    FILE *fp = NULL;
    bool is_error = false;

	if ((fp = fopen(file_name, "r")) == NULL) {
		printf("\topp_load_from_file - Unable to open file %s\n", file_name);
		exit(-1);
	}
	if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2) {
		printf("\topp_load_from_file - error reading file data from %s\n", file_name);
		exit(-1);
	}
    if (fsize < set_size || fdim != dim) {
		printf("\topp_load_from_file - dim and/or set_size issue in file %s\n", file_name);
		exit(-1);        
    }

    void* data_ptr = (void *)opp_host_malloc((size_t)(set_size * dim * size));

    if (strncmp("double", type, 6) == 0) {
        double* data = (double*)data_ptr;

        for (int n = 0; n < set_size; n++) {
            switch (dim) {
            case 1:
                if (fscanf(fp, " %lf\n", &data[n * dim + 0]) != 1) 
                    is_error = true;
                break;
            case 2:
                if (fscanf(fp, " %lf %lf\n", &data[n * dim + 0], &data[n * dim + 1]) != 2) 
                    is_error = true;
                break;
            case 3:
                if (fscanf(fp, " %lf %lf %lf\n", 
                        &data[n * dim + 0], &data[n * dim + 1], &data[n * dim + 2]) != 3) 
                    is_error = true;
                break;
            case 4:
                if (fscanf(fp, " %lf %lf %lf %lf\n", 
                        &data[n * dim + 0], &data[n * dim + 1], &data[n * dim + 2], &data[n * dim + 3]) != 4) 
                    is_error = true;
                break;    
            case 16:
                if (fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                        &data[n * dim + 0], &data[n * dim + 1], &data[n * dim + 2], &data[n * dim + 3], 
                        &data[n * dim + 4], &data[n * dim + 5], &data[n * dim + 6], &data[n * dim + 7],
                        &data[n * dim + 8], &data[n * dim + 9], &data[n * dim + 10], &data[n * dim + 11], 
                        &data[n * dim + 12], &data[n * dim + 13], &data[n * dim + 14], &data[n * dim + 15]) != 16) 
                    is_error = true;
                break;
            default: is_error = true;
            }
                    
            if (is_error) {
                opp_printf("opp_load_from_file", "Error reading Double from %s at index %d", file_name, n);
                opp_host_free(data);
                exit(-1);
            }
        }
    } 
    else if (strncmp("int", type, 3) == 0) {
        int* data = (int*)data_ptr;

        for (int n = 0; n < set_size; n++) {
            switch (dim) {
            case 1:
                if (fscanf(fp, " %d\n", &data[n * dim + 0]) != 1) 
                    is_error = true;
                break;
            case 2:
                if (fscanf(fp, " %d %d\n", &data[n * dim + 0], &data[n * dim + 1]) != 2) 
                    is_error = true;
                break;
            case 3:
                if (fscanf(fp, " %d %d %d\n", 
                    &data[n * dim + 0], &data[n * dim + 1], &data[n * dim + 2]) != 3) 
                        is_error = true;
                break;
            case 4:
                if (fscanf(fp, " %d %d %d %d\n", 
                    &data[n * dim + 0], &data[n * dim + 1], &data[n * dim + 2], &data[n * dim + 3]) != 4) 
                        is_error = true;
                break;    
            default: is_error = true;
            }
                    
            if (is_error) {
                opp_printf("opp_load_from_file", "Error reading Int from %s at index %d", file_name, n);
                opp_host_free(data);
                exit(-1);
            }
        }
    } 
    else {
        opp_printf("opp_load_from_file", "Error Unsupported type for loading %s", type);
        opp_host_free(data_ptr);
        exit(-2);
    }

    return data_ptr;
}
