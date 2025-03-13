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
#include <opp_mpi.h>

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
    MPI_CHECK(MPI_Win_free(&win_structMeshToCellMapping));
    structMeshToCellMapping = nullptr;

    MPI_CHECK(MPI_Win_free(&win_structMeshToRankMapping));
    structMeshToRankMapping = nullptr;           
};

//*******************************************************************************
void CellMapper::createStructMeshMappingArrays() 
{ 
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    //*************************** // One per shared memory (node)
    auto createPerNodeSharedMemArrays = [&](int*& mapping, MPI_Win& win) {
   
        const MPI_Aint size = ((comm->rank_intra == 0) ? (globalGridSize * sizeof(int)) : 0);

        MPI_CHECK(MPI_Win_allocate_shared(size, sizeof(int), MPI_INFO_NULL, comm->comm_intra, 
                    (void *)&mapping, &win));

        MPI_Aint allocatedSize = 0;
        int disp = 0;
        MPI_CHECK(MPI_Win_shared_query(win, 0, &allocatedSize, &disp, (void *)&mapping));

        if (globalGridSize * sizeof(int) != (size_t)allocatedSize) {
            std::stringstream ss;
            ss << "Pointer to incorrect size in Mapping; globalGridSize:" << globalGridSize << 
                " allocatedSize:" << allocatedSize << " at " << __FILE__ << ":" << __LINE__; 
            opp_abort(ss.str());
        }
        if (disp != sizeof(int)) {
            std::stringstream ss;
            ss << "Invalid displacement in Mapping; disp:" << disp << " at " << __FILE__ << ":" << __LINE__; 
            opp_abort(ss.str());
        }

        if (comm->rank_intra == 0) {
            for (size_t i = 0; i < globalGridSize; i++)
                mapping[i] = MAX_CELL_INDEX;
        }

        MPI_CHECK(MPI_Win_fence(0, win));
    };

    createPerNodeSharedMemArrays(structMeshToCellMapping, win_structMeshToCellMapping);
    createPerNodeSharedMemArrays(structMeshToRankMapping, win_structMeshToRankMapping);

    waitBarrier();
}

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) 
{
    waitBarrier();

    if (comm->rank_intra == 0) {

        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToCellMapping));
        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToRankMapping));

        if (OPP_DBG) 
            opp_printf("CellMapper", "reduceInterNodeMappings Start size_inter %d", comm->size_inter);

        std::stringstream s_log("Send Counts: "), r_log("Recv Counts: ");

        int totalSent = 0, totalRecv = 0;
        std::vector<int> sendCounts(comm->size_inter);
        std::vector<int> sendDisplacements(comm->size_inter);
        std::vector<int> recvCounts(comm->size_inter);
        std::vector<int> recvDisplacements(comm->size_inter);

        for (int i = 0; i < comm->size_inter; ++i) {
            const int remainder = (int)(globalGridSize % comm->size_inter);
            sendCounts[i] = (globalGridSize / comm->size_inter) + (i < remainder ? 1 : 0);
            sendDisplacements[i] = totalSent;
            totalSent += sendCounts[i];

            if (OPP_DBG) 
                s_log << sendCounts[i] << "|Disp:" << sendDisplacements[i] << " ";
        }

        for (int i = 0; i < comm->size_inter; ++i) {
            recvCounts[i] = sendCounts[comm->rank_inter];
            recvDisplacements[i] = totalRecv;
            totalRecv += recvCounts[i];

            if (OPP_DBG) 
                r_log << recvCounts[i] << "|Disp:" << recvDisplacements[i] << " ";
        }

        if (OPP_DBG) {
            opp_printf("CellMapper::reduceInterNodeMappings", "SendCount: %d Disp: %s", 
                            totalSent, s_log.str().c_str());
            opp_printf("CellMapper::reduceInterNodeMappings", "RecvCount: %d Disp: %s", 
                            totalRecv, r_log.str().c_str());
        }

        std::vector<int> cellMappingsRecv(totalRecv, MAX_CELL_INDEX-1);
        std::vector<int> ranksRecv(totalRecv, MAX_CELL_INDEX-1);

        const int INTS_PER_BATCH = 10000;

        std::vector<MPI_Request> sendRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) {
            if (rank == comm->rank_inter) 
                continue;

            const int totalSendCount = sendCounts[rank];
            int sendCount = 0, alreadySent = 0;

            const int batches = (totalSendCount / INTS_PER_BATCH) + 
                                    ((totalSendCount % INTS_PER_BATCH == 0) ? 0 : 1);

            for (int i = 0; i < batches; i++) {
                
                sendCount = ((i != batches - 1) ? 
                                INTS_PER_BATCH : (totalSendCount - alreadySent));

                sendRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Isend(
                            &(structMeshToRankMapping[sendDisplacements[rank] + i * INTS_PER_BATCH]), 
                            sendCount, MPI_INT, rank, (10000 + i), 
                            comm->comm_inter, &sendRequests[sendRequests.size() - 1])); 

                sendRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Isend(
                            &(structMeshToCellMapping[sendDisplacements[rank] + i * INTS_PER_BATCH]), 
                            sendCount, MPI_INT, rank, (20000 + i), 
                            comm->comm_inter, &sendRequests[sendRequests.size() - 1])); 

                alreadySent += sendCount;
            }
        }

        std::vector<MPI_Request> recvRequests;
        for (int rank = 0; rank < comm->size_inter; rank++) {
            if (rank == comm->rank_inter) 
                continue;

            const int totalRecvCount = recvCounts[rank];
            int recvCount = 0, alreadyRecvd = 0;

            const int batches = (totalRecvCount / INTS_PER_BATCH) + 
                                    ((totalRecvCount % INTS_PER_BATCH == 0) ? 0 : 1);

            for (int i = 0; i < batches; i++) {

                recvCount = ((i != batches-1) ? 
                                INTS_PER_BATCH : (totalRecvCount - alreadyRecvd));

                // opp_printf("RECV", "batches %d|%d disp %d count %d from rank %d", i, batches,
                //     recvDisplacements[rank]+i*INTS_PER_BATCH, recvCount, inter_ranks[rank]);

                recvRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Irecv(
                            &(ranksRecv[recvDisplacements[rank] + i * INTS_PER_BATCH]), 
                            recvCount, MPI_INT, rank, (10000 + i), comm->comm_inter, 
                            &recvRequests[recvRequests.size() - 1]));

                recvRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Irecv(
                            &(cellMappingsRecv[recvDisplacements[rank] + i * INTS_PER_BATCH]), 
                            recvCount, MPI_INT, rank, (20000 + i), comm->comm_inter, 
                            &recvRequests[recvRequests.size() - 1]));

                alreadyRecvd += recvCount;
            }
        }

        // Copy own data
        for (int i = 0; i < recvCounts[comm->rank_inter]; i++) {
            
            const int recv_disp = recvDisplacements[comm->rank_inter];
            const int send_disp = sendDisplacements[comm->rank_inter];

            ranksRecv[recv_disp + i] = structMeshToRankMapping[send_disp + i];
            cellMappingsRecv[recv_disp + i] = structMeshToCellMapping[send_disp + i];
        }

        MPI_CHECK(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE));
        MPI_CHECK(MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE));
        MPI_CHECK(MPI_Barrier(comm->comm_inter));

        const size_t recvCount = recvCounts[0];
        for (size_t i = 0; i < recvCount; i++) { // reduce to get common mapping on all inter node ranks

            int cellIndex = MAX_CELL_INDEX;
            for (int r = 0; r < comm->size_inter; r++) {

                if (cellIndex > cellMappingsRecv[i + r * recvCount]) {

                    cellIndex = cellMappingsRecv[i + r * recvCount];
                    cellMappingsRecv[i] = cellIndex;
                    ranksRecv[i] = ranksRecv[i + r * recvCount];
                }
            }
        }

        MPI_CHECK(MPI_Barrier(comm->comm_inter));

        sendRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) {

            if (rank == comm->rank_inter) 
                continue;

            const int total_sent = recvCount;
            int sendCount = 0, alreadySent = 0;

            const int batches = (total_sent / INTS_PER_BATCH) + 
                                    ((total_sent % INTS_PER_BATCH == 0) ? 0 : 1);

            for (int i = 0; i < batches; i++) {
                
                sendCount = ((i != batches - 1) ? INTS_PER_BATCH : (total_sent - alreadySent));

                sendRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Isend(
                            &(cellMappingsRecv[i * INTS_PER_BATCH]), sendCount, MPI_INT, rank, 
                            (30000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1])); 

                sendRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Isend(
                            &(ranksRecv[i * INTS_PER_BATCH]), sendCount, MPI_INT, rank, 
                            (40000 + i), comm->comm_inter, &sendRequests[sendRequests.size() - 1])); 

                alreadySent += sendCount;
            }
        }

        // Since we are about to write to data mapped by windows, release the lock here
        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToCellMapping));
        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToRankMapping));

        recvRequests.clear();
        for (int rank = 0; rank < comm->size_inter; rank++) {

            if (rank == comm->rank_inter) 
                continue;

            const int total_recv = sendCounts[rank];
            int recvCount = 0, alreadyRecvd = 0;

            const int batches = (total_recv / INTS_PER_BATCH) + 
                                    ((total_recv % INTS_PER_BATCH == 0) ? 0 : 1);

            for (int i = 0; i < batches; i++) {

                recvCount = ((i != batches-1) ? INTS_PER_BATCH : (total_recv - alreadyRecvd));

                recvRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Irecv(
                            &(structMeshToCellMapping[sendDisplacements[rank] + i * INTS_PER_BATCH]), 
                            recvCount, MPI_INT, rank, (30000 + i), comm->comm_inter, 
                            &recvRequests[recvRequests.size() - 1]));

                recvRequests.emplace_back(MPI_Request());
                MPI_CHECK(MPI_Irecv(
                            &(structMeshToRankMapping[sendDisplacements[rank] + i * INTS_PER_BATCH]), 
                            recvCount, MPI_INT, rank, (40000 + i), comm->comm_inter, 
                            &recvRequests[recvRequests.size() - 1]));

                alreadyRecvd += recvCount;
            }
        }

        // Copy own data
        for (int i = 0; i < sendCounts[comm->rank_inter]; i++) {
            
            const int send_disp = sendDisplacements[comm->rank_inter];
            structMeshToRankMapping[send_disp + i] = ranksRecv[i];
            structMeshToCellMapping[send_disp + i] = cellMappingsRecv[i];
        }

        MPI_CHECK(MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUS_IGNORE));
        MPI_CHECK(MPI_Waitall(recvRequests.size(), recvRequests.data(), MPI_STATUS_IGNORE));
        MPI_CHECK(MPI_Barrier(comm->comm_inter));
        
        if (OPP_DBG) 
            opp_printf("CellMapper", "reduceInterNodeMappings END");
    }

    waitBarrier();
}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat);

    if (comm->rank_intra == 0) {    
        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToCellMapping));

        for (size_t i = 0; i < globalGridSize; i++) {         
            if (structMeshToCellMapping[i] != MAX_CELL_INDEX)
                structMeshToCellMapping[i] = (-1 * structMeshToCellMapping[i]);
        }

        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToCellMapping));
    }

    waitBarrier();

    convertToLocalMappings_seq(globalToLocalCellIndexMapper);

    MPI_CHECK(MPI_Barrier(OPP_MPI_WORLD));

    if (comm->rank_intra == 0) {
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToCellMapping, globalGridSize, 
                        MPI_INT, MPI_MAX, comm->comm_inter));
    }

    waitBarrier();

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::convertToLocalMappingsIncRank(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappingsIncRank Start");

    GlobalToLocalCellIndexMapper globalToLocalCellIndexMapper(global_cell_id_dat, false);

    if (comm->rank_intra == 0) {    
        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToCellMapping));

        for (size_t i = 0; i < globalGridSize; i++) {         
            if (structMeshToCellMapping[i] != MAX_CELL_INDEX) {
                structMeshToCellMapping[i] = (-1 * structMeshToCellMapping[i]);
            }
        }

        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToCellMapping));
    }

    waitBarrier();

    std::vector<int> tmp_structMeshToCellMapping(globalGridSize);
    std::copy(structMeshToCellMapping, structMeshToCellMapping + globalGridSize, tmp_structMeshToCellMapping.begin());

    for (size_t i = 0; i < globalGridSize; i++) {

        const int globalCID = (-1 * tmp_structMeshToCellMapping[i]);
        if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {               
            
            const int localCID = globalToLocalCellIndexMapper.map(globalCID);   
            if (localCID != MAX_CELL_INDEX) {
                enrichStructuredMesh(i, localCID, OPP_rank);
            }
        }
    }

    waitBarrier();

    if (comm->rank_intra == 0) {
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToCellMapping, globalGridSize, 
                        MPI_INT, MPI_MAX, comm->comm_inter));
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToRankMapping, globalGridSize, 
                        MPI_INT, MPI_MIN, comm->comm_inter));
    }

    waitBarrier();

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappingsIncRank END");
}

//*******************************************************************************
void CellMapper::generateStructuredMeshFromFile(opp_set set, const opp_dat c_gbl_id) {

    if (OPP_rank == 0)            
        opp_printf("OPP", "generateStructuredMeshFromFile START cells [%s] global grid dims %zu %zu %zu",
            set->name, globalGridDimsX, globalGridDimsY, globalGridDimsZ);

    createStructMeshMappingArrays();

    opp_profiler->start("Setup_Mover_s0");

    int set_size = 0;
    MPI_Allreduce(&(set->size), &set_size, 1, MPI_INT, MPI_SUM, OPP_MPI_WORLD);

    if (comm->rank_intra == 0) { // read only with the node-main rank   
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

    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance,
    //      however, not like in generating values, at this time we dont have structMeshToRankMapping enriched!
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
    convertToLocalMappingsIncRank(c_gbl_id);
    opp_profiler->end("Setup_Mover_s5");

    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile DONE");
}

//*******************************************************************************
// Use OPP_DH_DATA_DUMP=1 to dump the structured mesh
void CellMapper::generateStructuredMesh(opp_set set, const opp_dat c_gbl_id,
            const std::function<void(const opp_point&, int&)>& all_cell_checker) { 

    if (OPP_rank == 0)            
        opp_printf("OPP", "generateStructuredMesh START cells [%s] global grid dims %zu %zu %zu DUMP %s",
            set->name, globalGridDimsX, globalGridDimsY, globalGridDimsZ, OPP_dh_data_dump? "YES" : "NO");

    const int set_size_inc_halo = set->size + set->exec_size + set->nonexec_size;
    if (set_size_inc_halo <= 0) {
        opp_printf("OPP", "Error... set_size_inc_halo <= 0 for set %s", set->name);
        opp_abort("Error... APP set_size_inc_halo <= 0");
    }

    std::map<size_t, opp_point> removed_coords;
    const opp_point& min_glb_coords = boundingBox->getGlobalMin();
    const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT define

    createStructMeshMappingArrays();

    // Step 1 : Get centroids of the structured mesh cells and try to relate them to unstructured mesh indices
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 1 Start");
    opp_profiler->start("Setup_Mover_s1");

    for (size_t dz = localGridStart.z; dz < localGridEnd.z; dz++) {       
        double z = min_glb_coords.z + dz * gridSpacing;        
        for (size_t dy = localGridStart.y; dy < localGridEnd.y; dy++) {            
            double y = min_glb_coords.y + dy * gridSpacing;           
            for (size_t dx = localGridStart.x; dx < localGridEnd.x; dx++) {                
                double x = min_glb_coords.x + dx * gridSpacing;               
                
                size_t index = (dx + dy * globalGridDimsX + dz * globalGridDimsXY); 

                opp_point centroid = getCentroidOfBox(opp_point(x, y ,z));
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
                    enrichStructuredMesh(index, ((int*)c_gbl_id->data)[cid], OPP_rank);
                }
            }
        }
    }
    opp_profiler->end("Setup_Mover_s1");

    // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 2 Start");
    opp_profiler->start("Setup_Mover_s2");
    reduceInterNodeMappings(1);

    // The marked structured cells from this rank might be filled by another rank, so if already filled, 
    // no need to recalculate from current rank
    for (auto it = removed_coords.begin(); it != removed_coords.end(); ) {
        size_t removed_idx = it->first;
        if (structMeshToRankMapping[removed_idx] != MAX_CELL_INDEX) {
            it = removed_coords.erase(it); // This structured index is already written by another rank
            // opp_printf("OPP", "index %zu already in %d", this->structMeshToRankMapping[removed_idx], removed_idx);
        } 
        else {
            ++it;
        } 
    }
    waitBarrier();    
    opp_profiler->end("Setup_Mover_s2");

    // Step 3 : Iterate over NEED_REMOVE points, Check whether atleast one vertex of the structured mesh is within 
    //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 3 Start");
    opp_profiler->start("Setup_Mover_s3");
    for (auto& p : removed_coords) {

        const size_t index = p.first;
        double &x = p.second.x, &y = p.second.y, &z = p.second.z;
        
        const double gs = gridSpacing;
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
        lockWindows();
        const int avail_gbl_cid = structMeshToCellMapping[index]; 
        if ((most_suitable_gbl_cid != MAX_CELL_INDEX) && (most_suitable_gbl_cid < avail_gbl_cid) && 
                    (most_suitable_cid < set->size)) {
            enrichStructuredMesh(index, most_suitable_gbl_cid, OPP_rank);   
        }
        unlockWindows();
    }
    opp_profiler->end("Setup_Mover_s3");

    // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 4 Start");
    opp_profiler->start("Setup_Mover_s4");
    reduceInterNodeMappings(2);
    opp_profiler->end("Setup_Mover_s4");

    // Step Add : Dump the structured mesh to a file, if requested 
    if (OPP_dh_data_dump) {
        int set_size = 0;
        MPI_Reduce(&(set->size), &set_size, 1, MPI_INT, MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);

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
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
    convertToLocalMappings(c_gbl_id);
    opp_profiler->end("Setup_Mover_s5");

    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh DONE");
}
