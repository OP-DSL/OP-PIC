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
#include <opp_cuda.h>

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

    opp_mem::copy_host_to_dev<OPP_REAL>(oneOverGridSpacing_d, &oneOverGridSpacing, 1, false, true, 1);
    
    opp_mem::copy_host_to_dev<OPP_REAL>(minGlbCoordinate_d, (OPP_REAL*)&(minGlbCoordinate), 
                                        dim, false, true, dim);
    
    const size_t temp[4] = { globalGridDimsX, globalGridDimsY, globalGridDimsZ, globalGridDimsXY };
    opp_mem::copy_host_to_dev<size_t>(globalGridDims_d, temp, 4, false, true, 4);

    opp_mem::copy_host_to_dev<size_t>(globalGridSize_d, &globalGridSize, 1, false, true, 1);
}

//*******************************************************************************
CellMapper::~CellMapper() 
{ 
#ifdef USE_MPI
    MPI_CHECK(MPI_Win_free(&win_structMeshToCellMapping));
    structMeshToCellMapping = nullptr;

    MPI_CHECK(MPI_Win_free(&win_structMeshToRankMapping));
    structMeshToRankMapping = nullptr;  

    opp_mem::dev_free(structMeshToRankMapping_d);         
#else
    delete[] structMeshToCellMapping;
    structMeshToCellMapping = nullptr;
#endif

    opp_mem::dev_free(structMeshToCellMapping_d);
};

//*******************************************************************************
void CellMapper::createStructMeshMappingArrays() 
{
#ifdef USE_MPI 
    MPI_CHECK(MPI_Barrier(OPP_MPI_WORLD));

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

    MPI_CHECK(MPI_Barrier(OPP_MPI_WORLD));

    structMeshToRankMapping_d = opp_mem::dev_malloc<OPP_INT>(globalGridSize);
#else
    structMeshToCellMapping = new int[globalGridSize];
    for (size_t i = 0; i < globalGridSize; i++)
        structMeshToCellMapping[i] = MAX_CELL_INDEX;
#endif

    structMeshToCellMapping_d = opp_mem::dev_malloc<OPP_INT>(globalGridSize);
}

//*******************************************************************************
void CellMapper::reduceInterNodeMappings(int callID) 
{
#ifdef USE_MPI
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
#endif
}

//*******************************************************************************
void CellMapper::convertToLocalMappings(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings Start");

#ifdef USE_MPI
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

    waitBarrier();

    if (comm->rank_intra == 0) {
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToCellMapping, globalGridSize, 
                        MPI_INT, MPI_MAX, comm->comm_inter));
    }

    waitBarrier();
#endif
    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappings END");
}

//*******************************************************************************
void CellMapper::hostToDeviceTransferMappings() {

    if (OPP_DBG) 
        opp_printf("CellMapper", "hostToDeviceTransferMappings Start");

    opp_mem::copy_host_to_dev<OPP_INT>(structMeshToCellMapping_d, 
                            structMeshToCellMapping, globalGridSize);

#ifdef USE_MPI
    opp_mem::copy_host_to_dev<OPP_INT>(structMeshToRankMapping_d, 
                            structMeshToRankMapping, globalGridSize);
#endif

    if (OPP_DBG) {
        std::vector<OPP_INT> tmp_cells(globalGridSize), tmp_ranks(globalGridSize);

        opp_mem::copy_dev_to_host<OPP_INT>(tmp_cells.data(), 
                            structMeshToCellMapping_d, globalGridSize);
        opp_mem::copy_dev_to_host<OPP_INT>(tmp_ranks.data(), 
                            structMeshToRankMapping_d, globalGridSize);

        // for (int r = 0; r < OPP_comm_size; r++) {
        //     if (r == OPP_rank) {
        //         printStructuredMesh("structMeshToCellMappingBC", tmp_cells.data(), 
        //                                 globalGridSize, false, globalGridDimsX);
        //         printStructuredMesh("structMeshToRankMappingBC", tmp_ranks.data(), 
        //                                 globalGridSize, false, globalGridDimsX);   
        //     }
        //     MPI_Barrier(OPP_MPI_WORLD);  
        // }
    }

    if (OPP_DBG) 
        opp_printf("CellMapper", "hostToDeviceTransferMappings END");
}

//*******************************************************************************
void CellMapper::generateStructuredMesh(opp_set set, const opp_dat c_gbl_id, 
            const std::function<void(const opp_point&, int&)>& all_cell_checker) { 

    if (OPP_rank == 0)            
        opp_printf("OPP", "generateStructuredMesh START cells [%s] global grid dims %zu %zu %zu",
            set->name, globalGridDimsX, globalGridDimsY, globalGridDimsZ);

    const int set_size_inc_halo = set->size + set->exec_size + set->nonexec_size;
    if (set_size_inc_halo <= 0) {
        opp_printf("OPP", "Error... set_size_inc_halo <= 0 for set %s", set->name);
        opp_abort("Error... APP set_size_inc_halo <= 0");
    }

    std::map<size_t, opp_point> removed_coords;
    const opp_point& min_glb_coords = boundingBox->getGlobalMin();
    const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT define

    createStructMeshMappingArrays();

#ifdef USE_OMP
    const int omp_nthreads = omp_get_max_threads();
    opp_printf("OPP", "generateStructuredMesh omp_nthreads %d", omp_nthreads);
#endif

    // Step 1 : Get centroids of the structured mesh cells and try to relate them to unstructured mesh indices
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 1 Start");
    opp_profiler->start("Setup_Mover_s1");

#ifdef USE_OMP
    #pragma omp parallel for
    for (int thr = 0; thr < omp_nthreads; thr++)
    {
        const size_t start  = ((localGridEnd.z - localGridStart.z) * thr) / omp_nthreads;
        const size_t end = ((localGridEnd.z - localGridStart.z) * (thr+1)) / omp_nthreads;
        
        opp_printf("OPP", "generateStructuredMesh Starting thr %d [%zu %zu]", thr, start, end);

        for (size_t dz = start; dz < end; dz++)
        {

#else   
    for (size_t dz = localGridStart.z; dz < localGridEnd.z; dz++) {     
#endif  
        double z = min_glb_coords.z + dz * gridSpacing;        
        for (size_t dy = localGridStart.y; dy < localGridEnd.y; dy++) {            
            double y = min_glb_coords.y + dy * gridSpacing;
            for (size_t dx = localGridStart.x; dx < localGridEnd.x; dx++) {                
                double x = min_glb_coords.x + dx * gridSpacing;               
                
                size_t index = (dx + dy * globalGridDimsX + dz * globalGridDimsXY); 

                const opp_point centroid = getCentroidOfBox(opp_point(x, y ,z));
                int cid = MAX_CELL_INDEX;

                all_cell_checker(centroid, cid); // Find in which cell this centroid lies

                if (cid == MAX_CELL_INDEX) {
#ifdef USE_OMP
                    #pragma omp critical
                    {
                        removed_coords.insert(std::make_pair(index, opp_point(x, y ,z)));
                    }
#else
                    removed_coords.insert(std::make_pair(index, opp_point(x, y ,z)));
#endif
                }
                else if (cid < set->size) { // write only if the structured cell belong to the current MPI rank                    
                    enrichStructuredMesh(index, ((int*)c_gbl_id->data)[cid], OPP_rank);
                }
            }
        }
#ifdef USE_OMP
    }
#endif
    }
    opp_profiler->end("Setup_Mover_s1");

    // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 2 Start");
    opp_profiler->start("Setup_Mover_s2");
#ifdef USE_MPI
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
#endif
    opp_profiler->end("Setup_Mover_s2");

    // Step 3 : Iterate over NEED_REMOVE points, Check whether atleast one vertex of the structured mesh is within 
    //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 3 Start");
    opp_profiler->start("Setup_Mover_s3");
#ifdef USE_OMP 
    std::vector<size_t> removed_coords_keys;
    removed_coords_keys.reserve(removed_coords.size());
    for (const auto& pair : removed_coords)
        removed_coords_keys.push_back(pair.first);

    std::vector<std::vector<std::pair<int, int>>> tmp_add_per_thr;
    tmp_add_per_thr.resize(omp_nthreads);

    #pragma omp parallel for
    for (int thr = 0; thr < omp_nthreads; thr++)
    {
        const size_t start  = (removed_coords_keys.size() * thr) / omp_nthreads;
        const size_t end = (removed_coords_keys.size() * (thr+1)) / omp_nthreads;
      
        for (size_t i = start; i < end; i++)
        {

        const size_t index = removed_coords_keys[i];
        opp_point& p = removed_coords[index];
        double &x = p.x, &y = p.y, &z = p.z;
#else
    for (auto& p : removed_coords) {

        const size_t index = p.first;
        double &x = p.second.x, &y = p.second.y, &z = p.second.z;
#endif
        
        const double gs = gridSpacing;
        int most_suitable_cid = MAX_CELL_INDEX, most_suitable_gbl_cid = MAX_CELL_INDEX;

        std::array<opp_point,4> vertices = {
            opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z)),
        };

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
#ifndef USE_OMP
        lockWindows();
#endif
        const int avail_gbl_cid = structMeshToCellMapping[index]; 
        if ((most_suitable_gbl_cid != MAX_CELL_INDEX) && (most_suitable_gbl_cid < avail_gbl_cid) && 
                    (most_suitable_cid < set->size)) {
#ifdef USE_OMP
            tmp_add_per_thr[thr].push_back(std::make_pair(index, most_suitable_gbl_cid));
#else
            enrichStructuredMesh(index, most_suitable_gbl_cid, OPP_rank);
#endif     
        }
#ifndef USE_OMP
        unlockWindows();
#else
    }
#endif
    }
    opp_profiler->end("Setup_Mover_s3");

    // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 4 Start");
    opp_profiler->start("Setup_Mover_s4");
    reduceInterNodeMappings(2);
    opp_profiler->end("Setup_Mover_s4");

    // Step Add : Dump the structured mesh to a file, if requested 
    if (OPP_dh_data_dump) {
        opp_printf("OPP", "generateStructuredMesh: Not Dumping DH file, Use a CPU backend to Dump");
    }

    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
    convertToLocalMappings(c_gbl_id);
    opp_profiler->end("Setup_Mover_s5");

    // Step 6 : Copy structured mesh mappings from host to device
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 6 Start");
    opp_profiler->start("Setup_Mover_s6");
    hostToDeviceTransferMappings();
    opp_profiler->end("Setup_Mover_s6");

    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh DONE");
}


//*******************************************************************************
void CellMapper::convertToLocalMappingsIncRank(const opp_dat global_cell_id_dat) {

    if (OPP_DBG) 
        opp_printf("CellMapper", "convertToLocalMappingsIncRank Start");

#ifdef USE_MPI
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

    for (size_t i = 0; i < globalGridSize; i++) {

        const int globalCID = (-1 * structMeshToCellMapping[i]);
        if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {               
            
            const int localCID = globalToLocalCellIndexMapper.map(globalCID);   
            if (localCID != MAX_CELL_INDEX) {
                enrichStructuredMesh(i, localCID, OPP_rank);
            }
        }
    }

    MPI_CHECK(MPI_Barrier(OPP_MPI_WORLD));

    if (comm->rank_intra == 0) {
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToCellMapping, globalGridSize, 
                        MPI_INT, MPI_MAX, comm->comm_inter));
        MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE, structMeshToRankMapping, globalGridSize, 
            MPI_INT, MPI_MIN, comm->comm_inter));
    }

    waitBarrier();
#endif

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

#ifdef USE_MPI
    int set_size = 0;
    MPI_Allreduce(&(set->size), &set_size, 1, MPI_INT, MPI_SUM, OPP_MPI_WORLD);

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

    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance,
    //      however, not like in generating values, at this time we dont have structMeshToRankMapping enriched!
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
    convertToLocalMappingsIncRank(c_gbl_id);
    opp_profiler->end("Setup_Mover_s5");

    // Step 6 : Copy structured mesh mappings from host to device
    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMesh Step 6 Start");
    opp_profiler->start("Setup_Mover_s6");
    hostToDeviceTransferMappings();
    opp_profiler->end("Setup_Mover_s6");

    if (OPP_rank == 0) opp_printf("OPP", "generateStructuredMeshFromFile DONE");
}

//*******************************************************************************
void opp_init_dh_device(opp_set set)
{
    if (OPP_DBG) opp_printf("OPP", "opp_init_dh_device START");

    // const OPP_INT buffer_alloc_size = (set->set_capacity / 4);
    const OPP_INT buffer_alloc_size = set->set_capacity;
    if (dh_indices_h.capacity < buffer_alloc_size) {
        
        if (OPP_DBG) opp_printf("OPP", "opp_init_dh_device Init arrays %d", buffer_alloc_size);

        dh_indices_h.capacity = buffer_alloc_size;
        dh_indices_d.capacity = buffer_alloc_size;

        opp_mem::host_free(dh_indices_h.part_indices); 
        opp_mem::host_free(dh_indices_h.cell_indices); 
        opp_mem::host_free(dh_indices_h.rank_indices); 

        opp_mem::dev_free(dh_indices_d.part_indices); 
        opp_mem::dev_free(dh_indices_d.cell_indices); 
        opp_mem::dev_free(dh_indices_d.rank_indices); 

        dh_indices_h.part_indices = opp_mem::host_malloc<OPP_INT>(dh_indices_h.capacity);
        dh_indices_h.cell_indices = opp_mem::host_malloc<OPP_INT>(dh_indices_h.capacity);
        dh_indices_h.rank_indices = opp_mem::host_malloc<OPP_INT>(dh_indices_h.capacity);

        dh_indices_d.part_indices = opp_mem::dev_malloc<OPP_INT>(dh_indices_d.capacity);
        dh_indices_d.cell_indices = opp_mem::dev_malloc<OPP_INT>(dh_indices_d.capacity);
        dh_indices_d.rank_indices = opp_mem::dev_malloc<OPP_INT>(dh_indices_d.capacity);

        if (dh_indices_h.move_count == nullptr)
            dh_indices_h.move_count = opp_mem::host_malloc<OPP_INT>(1);
        
        if (dh_indices_d.move_count == nullptr)
            dh_indices_d.move_count = opp_mem::dev_malloc<OPP_INT>(1);  
        
        if (OPP_DBG) opp_printf("OPP", "opp_init_dh_device capacities %d %d", 
                        dh_indices_h.capacity, dh_indices_d.capacity);    
    }

    *(dh_indices_h.move_count) = 0;
    opp_mem::copy_host_to_dev<OPP_INT>(dh_indices_d.move_count, dh_indices_h.move_count, 1);

    if (OPP_DBG) opp_printf("OPP", "opp_init_dh_device END");
}

#ifdef USE_MPI 
//*******************************************************************************
// gathers all global move information into the global mover for communication
void opp_gather_dh_move_indices(opp_set set)
{
    opp_profiler->start("MvDH_Gather");
    opp_mem::copy_dev_to_host<OPP_INT>(dh_indices_h.move_count, dh_indices_d.move_count, 1);

    if (OPP_DBG) 
        opp_printf("OPP", "opp_gather_dh_move_indices move_count %d", *(dh_indices_h.move_count));

    opp_mem::copy_dev_to_host<OPP_INT>(dh_indices_h.move_count, dh_indices_d.move_count, 1);

    opp_mem::copy_dev_to_host<OPP_INT>(dh_indices_h.part_indices, 
                                    dh_indices_d.part_indices, *(dh_indices_h.move_count));
    opp_mem::copy_dev_to_host<OPP_INT>(dh_indices_h.cell_indices, 
                                    dh_indices_d.cell_indices, *(dh_indices_h.move_count));
    opp_mem::copy_dev_to_host<OPP_INT>(dh_indices_h.rank_indices, 
                                    dh_indices_d.rank_indices, *(dh_indices_h.move_count));

// for (int r = 0; r < OPP_comm_size; r++) {
//     if (r == OPP_rank) {
//         std::stringstream ss1;
//         for (size_t i = 0; i < *(dh_indices_h.move_count); i++) {
//             ss1 << dh_indices_h.part_indices[i] << ",";
//             if ((i + 1) % 100 == 0) 
//                 ss1 << "\n";
//         }
//         printf("\n[dh_indices_h.part_indices %d on RANK - %d]\n%s\n", *(dh_indices_h.move_count), OPP_rank, ss1.str().c_str());
//         // std::stringstream ss2;
//         // for (size_t i = 0; i < *(dh_indices_h.move_count); i++) {
//         //     ss2 << dh_indices_h.rank_indices[i] << ",";
//         //     if ((i + 1) % 100 == 0) 
//         //         ss2 << "\n";
//         // }   
//         // printf("[dh_indices_h.rank_indices %d on RANK - %d]\n%s\n", *(dh_indices_h.move_count), OPP_rank, ss2.str().c_str());
//         std::stringstream ss3;
//         for (size_t i = 0; i < *(dh_indices_h.move_count); i++) {
//             ss3 << dh_indices_h.cell_indices[i] << ",";
//             if ((i + 1) % 100 == 0) 
//                 ss3 << "\n";
//         }   
//         printf("[dh_indices_h.cell_indices %d on RANK - %d]\n%s\n\n", *(dh_indices_h.move_count), OPP_rank, ss3.str().c_str());
//     }
//     MPI_Barrier(OPP_MPI_WORLD);
// }

    for (OPP_INT i = 0; i < *(dh_indices_h.move_count); i++) {
        globalMover->markParticleToMove(set, dh_indices_h.part_indices[i],
                            dh_indices_h.rank_indices[i], dh_indices_h.cell_indices[i]);
    }

    opp_profiler->end("MvDH_Gather");
    if (OPP_DBG) opp_printf("OPP", "opp_gather_dh_move_indices DONE - dh move count %d", 
                    *(dh_indices_h.move_count));
}

//*******************************************************************************
dh_particle_packer_gpu::dh_particle_packer_gpu(std::map<int, std::vector<OPP_INT>>& local_part_indices, 
                                            std::map<int, std::vector<OPP_INT>>& foreign_cell_indices) 
        : dh_particle_packer(local_part_indices, foreign_cell_indices)
{
    // opp_printf("dh_particle_packer_gpu::dh_particle_packer_gpu", "No implementation");
}

//*******************************************************************************
dh_particle_packer_gpu::~dh_particle_packer_gpu()
{
    // opp_printf("dh_particle_packer_gpu::~dh_particle_packer_gpu", "No implementation");
}

//*******************************************************************************
void dh_particle_packer_gpu::pack(opp_set set)
{
    if (OPP_DBG) 
        opp_printf("dh_particle_packer_gpu", "pack set [%s]", set->name);

    opp_profiler->start("MvDH_Pack");

    std::map<int, std::vector<char>>& buffers_of_set = this->buffers[set->index];
    thrust::device_vector<OPP_INT>& temp_dv = *(set->mesh_relation_dat->thrust_int_sort);

    for (const auto& per_rank_parts : local_part_ids) {

        const int send_rank = per_rank_parts.first;
        const std::vector<int>& part_ids_vec = per_rank_parts.second;
        const size_t bytes_per_rank = (size_t)set->particle_size * part_ids_vec.size();
        
        if (OPP_DBG) 
            opp_printf("dh_particle_packer_gpu", "pack send_rank %d - count %zu", 
                        send_rank, part_ids_vec.size());

        std::vector<char>& send_rank_buffer = buffers_of_set[send_rank];
        send_rank_buffer.resize(bytes_per_rank);

        const int copy_count = (int)part_ids_vec.size();
        
        OPP_INT* temp_dp = opp_get_dev_raw_ptr<OPP_INT>(temp_dv);
        opp_mem::copy_host_to_dev<OPP_INT>(temp_dp, part_ids_vec.data(), copy_count);

        size_t disp = 0;
        for (auto& dat : *(set->particle_dats)) {

            const size_t bytes_to_copy = (dat->size * part_ids_vec.size());

            if (dat->is_cell_index) {  
                memcpy(&(send_rank_buffer[disp]), foreign_cell_ids[send_rank].data(), bytes_to_copy);
            }
            else if (strcmp(dat->type, "double") == 0) {
                    
                copy_according_to_index<OPP_REAL>(dat->thrust_real, dat->thrust_real_sort, 
                    temp_dv, dat->set->set_capacity, copy_count, copy_count, dat->dim);

                opp_mem::copy_dev_to_host<char>(&(send_rank_buffer[disp]), 
                    (char*)opp_get_dev_raw_ptr<OPP_REAL>(*(dat->thrust_real_sort)), bytes_to_copy);
            }
            else if (strcmp(dat->type, "int") == 0) {
                
                copy_according_to_index<OPP_INT>(dat->thrust_int, dat->thrust_int_sort, 
                    temp_dv, dat->set->set_capacity, copy_count, copy_count, dat->dim);
                
                opp_mem::copy_dev_to_host<char>(&(send_rank_buffer[disp]), 
                    (char*)opp_get_dev_raw_ptr<OPP_INT>(*(dat->thrust_int_sort)), bytes_to_copy);
            }

            disp += bytes_to_copy; 
        }

        if (OPP_DBG)
            opp_printf("dh_particle_packer_gpu", "Packed %zu parts to send to rank %d, displacement %d", 
                send_rank_buffer.size(), send_rank, disp);
    }

    opp_profiler->end("MvDH_Pack");
}

//*******************************************************************************
void dh_particle_packer_gpu::unpack(opp_set set, const std::map<int, std::vector<char>>& part_recv_buffers,
                                    int64_t total_recv_count)
{
    if (OPP_DBG) 
        opp_printf("dh_particle_packer_gpu", "unpack - total_recv_count %lld set->size %d set->diff %d", 
            total_recv_count, set->size, set->diff);

    opp_profiler->start("MvDH_Unpack");
 
    if (total_recv_count > 0) {

        opp_increase_particle_count(set, (int)total_recv_count);

        std::vector<opp_dat>& particle_dats = *(set->particle_dats);

        // create a continuous memory in host to copy to device
        std::vector<std::vector<char>> temp_data_vec;
        for (const auto& dat : particle_dats) {
            temp_data_vec.push_back(std::vector<char>(dat->size * total_recv_count));
        }

        size_t current_recv_count = 0;
        for (const auto& a : part_recv_buffers) {

            const std::vector<char>& buffer = a.second;
            const size_t recv_count = (buffer.size() / set->particle_size);
            int64_t displacement = 0;

            for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++) {
            
                const opp_dat& dat = particle_dats[dat_idx];
                std::vector<char>& temp_data = temp_data_vec[dat_idx];

                const int element_size = dat->size / dat->dim;
                const size_t bytes_to_copy_per_dim = (element_size * recv_count);

                for (int i = 0; i < dat->dim; i++) {
                    memcpy(&(temp_data[element_size * (i * total_recv_count + current_recv_count)]), 
                        &(buffer[displacement]), bytes_to_copy_per_dim);

                    displacement += bytes_to_copy_per_dim; 
                }                
            }

            current_recv_count += recv_count;
        }

        const size_t recv_part_start_idx = (size_t)(set->size - set->diff);
        
        // copy to device
        for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++) {
        
            opp_dat dat = particle_dats[dat_idx];
            const std::vector<char>& temp_data = temp_data_vec[dat_idx];

            const size_t bytes_to_copy_per_dim = total_recv_count * dat->size / dat->dim;
            const int element_size = dat->size / dat->dim;

            for (size_t d = 0; d < (size_t)dat->dim; d++) {
                 
                const size_t data_d_offset = (recv_part_start_idx + d * set->set_capacity) * element_size;
                const size_t data_h_offset = d * total_recv_count * element_size;

                char* data_d = dat->data_d + data_d_offset;
                opp_mem::copy_host_to_dev<char>(data_d, &(temp_data[data_h_offset]), 
                                                bytes_to_copy_per_dim);   
            }
        }
    }

    opp_profiler->end("MvDH_Unpack");

    if (OPP_DBG) 
        opp_printf("dh_particle_packer_gpu", "Unpack END");
}
#endif