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

#include "opp_direct_hop_core.h"
#include <opp_seq.h>

//*******************************************************************************
using namespace opp;

std::shared_ptr<BoundingBox> boundingBox;
std::shared_ptr<CellMapper> cellMapper;
std::shared_ptr<Comm> comm;
std::unique_ptr<GlobalParticleMover> globalMover;
bool useGlobalMove = false;

//*******************************************************************************
CellMapper::CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, 
                        const std::shared_ptr<Comm> comm) 
    : boundingBox(boundingBox), gridSpacing(gridSpacing), oneOverGridSpacing(1.0 / gridSpacing), 
        minGlbCoordinate(boundingBox->getGlobalMin()), comm(comm), dim(boundingBox->getDim())
{
    init_host(gridSpacing);
}

//*******************************************************************************
CellMapper::~CellMapper() 
{ 
    delete[] structMeshToCellMapping;
    structMeshToCellMapping = nullptr;
};

//*******************************************************************************
void CellMapper::createStructMeshMappingArrays() 
{
    structMeshToCellMapping = new int[globalGridSize];
    for (size_t i = 0; i < globalGridSize; i++)
        structMeshToCellMapping[i] = MAX_CELL_INDEX;
}

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) 
{

}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat);

    convertToLocalMappings_seq(globalToLocalCellIndexMapper);

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::generateStructuredMeshFromFile(opp_set set, const opp_dat c_gbl_id) {

    if (OPP_rank == 0)            
        opp_printf("OPP", "generateStructuredMeshFromFile START cells [%s] global grid dims %zu %zu %zu",
            set->name, globalGridDimsX, globalGridDimsY, globalGridDimsZ);

    createStructMeshMappingArrays();

    opp_profiler->start("Setup_Mover_s0");

#ifdef USE_MPI
    int set_size = 0;
    MPI_Reduce(&(set->size), &set_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (comm->rank_intra == 0) // read only with the node-main rank
#else
    int set_size = set->size;
#endif
    {
        std::stringstream s;
        if (opp_params->get<OPP_STRING>("opp_dh_mesh_folder", false) != "NOT FOUND") {
            s << opp_params->get<OPP_STRING>("opp_dh_mesh_folder", false);
        }
        s << str(set_size, "dh_c%d");
        s << str(gridSpacing, "_gs%2.4lE");
        s << str(boundingBox->domain_expansion.x, "_ex%2.4lE");
        s << str(boundingBox->domain_expansion.y, "_%2.4lE");
        s << str(boundingBox->domain_expansion.z, "_%2.4lE.bin");

        opp_decompress_read(s.str(), globalGridSize * sizeof(int), structMeshToCellMapping);
    }
    opp_profiler->end("Setup_Mover_s0");

    // // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance,
    // //      however, not like in generating values, at this time we dont have structMeshToRankMapping enriched!
    // if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile Step 5 Start");
    // opp_profiler->start("Setup_Mover_s5");
    // convertToLocalMappingsIncRank(c_gbl_id);
    // opp_profiler->end("Setup_Mover_s5");

    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile DONE");
}

//*******************************************************************************
// Use OPP_DH_DATA_DUMP=1 to dump the structured mesh
void CellMapper::generateStructuredMesh(opp_set set, const opp_dat c_gbl_id, 
            const std::function<void(const opp_point&, int&)>& all_cell_checker) { 

    if (OPP_rank == 0)            
        opp_printf("APP", "generateStructuredMesh START cells [%s] global grid dims %zu %zu %zu",
            set->name, cellMapper->globalGridDimsX, cellMapper->globalGridDimsY, cellMapper->globalGridDimsZ);

    const int set_size_inc_halo = set->size + set->exec_size + set->nonexec_size;
    if (set_size_inc_halo <= 0) {
        opp_printf("APP", "Error... set_size_inc_halo <= 0 for set %s", set->name);
        opp_abort("Error... APP set_size_inc_halo <= 0");
    }

    std::map<size_t, opp_point> removed_coords;
    const opp_point& min_glb_coords = boundingBox->getGlobalMin();
    const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT define

    cellMapper->createStructMeshMappingArrays();

    // Step 1 : Get centroids of the structured mesh cells and try to relate them to unstructured mesh indices
    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh Step 1 Start");
    opp_profiler->start("Setup_Mover_s1");
    double x = 0.0, y = 0.0, z = 0.0;
    
    for (size_t dz = cellMapper->localGridStart.z; dz < cellMapper->localGridEnd.z; dz++) {       
        z = min_glb_coords.z + dz * cellMapper->gridSpacing;        
        for (size_t dy = cellMapper->localGridStart.y; dy < cellMapper->localGridEnd.y; dy++) {            
            y = min_glb_coords.y + dy * cellMapper->gridSpacing;           
            for (size_t dx = cellMapper->localGridStart.x; dx < cellMapper->localGridEnd.x; dx++) {                
                x = min_glb_coords.x + dx * cellMapper->gridSpacing;               
                
                size_t index = (dx + dy * cellMapper->globalGridDimsX + dz * cellMapper->globalGridDimsXY); 

                opp_point centroid = cellMapper->getCentroidOfBox(opp_point(x, y ,z));
                int cid = MAX_CELL_INDEX;

                all_cell_checker(centroid, cid); // Find in which cell this centroid lies

                if (cid == MAX_CELL_INDEX) { // Change the centroid slightly to avoid edge cases
                    centroid.x += gridSpacing / 100.0;
                    all_cell_checker(centroid, cid);
                }
                if (cid == MAX_CELL_INDEX) { // Change the centroid slightly to avoid edge cases
                    centroid.y += gridSpacing / 100.0;
                    all_cell_checker(centroid, cid);
                }
                if (cid == MAX_CELL_INDEX) { // Change the centroid slightly to avoid edge cases
                    centroid.z += gridSpacing / 100.0;
                    all_cell_checker(centroid, cid);
                }

                if (cid == MAX_CELL_INDEX) {
                    removed_coords.insert(std::make_pair(index, opp_point(x, y ,z)));
                }
                else if (cid < set->size) { // write only if the structured cell belong to the current MPI rank                    
                    cellMapper->enrichStructuredMesh(index, ((int*)c_gbl_id->data)[cid], OPP_rank);
                }
            }
        }
    }
    opp_profiler->end("Setup_Mover_s1");

    // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh Step 2 Start");
    opp_profiler->start("Setup_Mover_s2");
#ifdef USE_MPI
    cellMapper->reduceInterNodeMappings(1);

    // The marked structured cells from this rank might be filled by another rank, so if already filled, 
    // no need to recalculate from current rank
    for (auto it = removed_coords.begin(); it != removed_coords.end(); ) {
        size_t removed_idx = it->first;
        if (cellMapper->structMeshToRankMapping[removed_idx] != MAX_CELL_INDEX) {
            it = removed_coords.erase(it); // This structured index is already written by another rank
            // opp_printf("APP", "index %zu already in %d", this->structMeshToRankMapping[removed_idx], removed_idx);
        } 
        else {
            ++it;
        } 
    }

    cellMapper->waitBarrier();    
#endif
    opp_profiler->end("Setup_Mover_s2");

    // Step 3 : Iterate over NEED_REMOVE points, Check whether atleast one vertex of the structured mesh is within 
    //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh Step 3 Start");
    opp_profiler->start("Setup_Mover_s3");
    for (auto& p : removed_coords) {

        const size_t index = p.first;
        double &x = p.second.x, &y = p.second.y, &z = p.second.z;
        
        const double gs = cellMapper->gridSpacing;
        int most_suitable_cid = MAX_CELL_INDEX, most_suitable_gbl_cid = MAX_CELL_INDEX;

        std::vector<opp_point> vertices;
        vertices.push_back(opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z)));
        vertices.push_back(opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z)));
        vertices.push_back(opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z)));
        vertices.push_back(opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z)));
        if (dim == 3) {
            vertices.push_back(opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z+gs)));
            vertices.push_back(opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z+gs)));
            vertices.push_back(opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z+gs)));
            vertices.push_back(opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z+gs)));
        }

        for (const auto& point : vertices) {
            int cid = MAX_CELL_INDEX;

            all_cell_checker(point, cid);

            if ((cid != MAX_CELL_INDEX) && (cid < set->size)) { 
                const int gbl_cid = ((OPP_INT*)c_gbl_id->data)[cid];
                if (most_suitable_gbl_cid > gbl_cid) {
                    most_suitable_gbl_cid = gbl_cid;
                    most_suitable_cid = cid;
                }
            }
        }    

        // Allow neighbours to write on-behalf of the current rank, to reduce issues
        cellMapper->lockWindows();
        const int avail_gbl_cid = cellMapper->structMeshToCellMapping[index]; 
        if ((most_suitable_gbl_cid != MAX_CELL_INDEX) && (most_suitable_gbl_cid < avail_gbl_cid) && 
                    (most_suitable_cid < set->size)) {        
            cellMapper->enrichStructuredMesh(index, most_suitable_gbl_cid, OPP_rank);      
        }
        cellMapper->unlockWindows();
    }
    opp_profiler->end("Setup_Mover_s3");

    // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh Step 4 Start");
    opp_profiler->start("Setup_Mover_s4");
#ifdef USE_MPI
    cellMapper->reduceInterNodeMappings(2);
#endif
    opp_profiler->end("Setup_Mover_s4");

    // Step Add : Dump the structured mesh to a file, if requested 
    if (OPP_dh_data_dump) {
        int set_size = 0;

#ifdef USE_MPI       
        MPI_Reduce(&(set->size), &set_size, 1, MPI_INT, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
#else
        set_size = set->size;
#endif

        if (OPP_rank == OPP_ROOT) {

            std::stringstream s;
            s << str(set_size, "dh_c%d");
            s << str(gridSpacing, "_gs%2.4lE");
            s << str(boundingBox->domain_expansion.x, "_ex%2.4lE");
            s << str(boundingBox->domain_expansion.y, "_%2.4lE");
            s << str(boundingBox->domain_expansion.z, "_%2.4lE.bin");

            opp_printf("OPP", "generateStructuredMesh Step Dumping File Start");

            opp_compress_write(s.str(), structMeshToCellMapping, 
                globalGridSize);
            
            if (OPP_DBG) {
                std::vector<int> decompressedData(globalGridSize);
                opp_decompress_read(s.str(), globalGridSize * sizeof(int), 
                                    decompressedData.data());
                
                for (size_t i = 0; i < globalGridSize; i++) {
                    if (decompressedData[i] != structMeshToCellMapping[i]) {
                        opp_printf("OPP", "Incorrect value from file at %d - file %d - system %d",
                            i, decompressedData[i], structMeshToCellMapping[i]);
                    }
                }
            }
        }
    }

    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance
    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
#ifdef USE_MPI
    cellMapper->convertToLocalMappings(c_gbl_id);
#endif
    opp_profiler->end("Setup_Mover_s5");

    if (OPP_rank == 0) opp_printf("APP", "generateStructuredMesh DONE");
}

