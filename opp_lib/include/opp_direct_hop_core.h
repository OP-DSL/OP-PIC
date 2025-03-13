
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

#include <opp_lib_core.h>

#define BOUND_OFFSET 1e-12

namespace opp {

//*******************************************************************************
class BoundingBox {

private:
    int dim = 0;

    opp_point boundingBox[2];       // index 0 is min, index 1 is max
    opp_point globalBoundingBox[2]; // index 0 is min, index 1 is max

public:
    opp_point domain_expansion;

    //***********************************
    BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate) : dim(dim) { OPP_RETURN_IF_INVALID_PROCESS;

        this->boundingBox[0] = minCoordinate;
        this->boundingBox[1] = maxCoordinate;

        generateGlobalBoundingBox(0);

        if (OPP_DBG)
            opp_printf("Local BoundingBox [provided]", 
                "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
                this->boundingBox[0].x, this->boundingBox[0].y, this->boundingBox[0].z, 
                this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);

    }

    //***********************************
    BoundingBox(const opp_dat node_pos_dat, int dim, const opp_point expansion) 
        : dim(dim), domain_expansion(expansion) { OPP_RETURN_IF_INVALID_PROCESS;

        if (dim != 3 && dim != 2) {
            opp_abort(std::string("For now, only 2D/3D BoundingBox is implemented"));
        }

        const double* node_pos_data = (const double*)node_pos_dat->data;
        const int node_count = node_pos_dat->set->size;
        // const int node_count = node_pos_dat->set->size + node_pos_dat->set->exec_size + 
        //                             node_pos_dat->set->nonexec_size;

        opp_point minCoordinate = opp_point(MAX_REAL, MAX_REAL, MAX_REAL);
        opp_point maxCoordinate = opp_point(MIN_REAL, MIN_REAL, MIN_REAL);

        if (dim == 2) {
            minCoordinate.z = 0;
            maxCoordinate.z = 0;
        }

        // make the bounding box even over the halo regions
        for (int i = 0; i < node_count; i++) {
            minCoordinate.x = std::min(node_pos_data[i * dim + 0], minCoordinate.x);
            maxCoordinate.x = std::max(node_pos_data[i * dim + 0], maxCoordinate.x);
            minCoordinate.y = std::min(node_pos_data[i * dim + 1], minCoordinate.y);
            maxCoordinate.y = std::max(node_pos_data[i * dim + 1], maxCoordinate.y);
            if (dim == 3) {
                minCoordinate.z = std::min(node_pos_data[i * dim + 2], minCoordinate.z);
                maxCoordinate.z = std::max(node_pos_data[i * dim + 2], maxCoordinate.z);
            }
        }

        // Why BOUND_OFFSET? to overcome overlapping == nodes
        if (node_count != 0) {
            minCoordinate.x -= (domain_expansion.x + BOUND_OFFSET);
            minCoordinate.y -= (domain_expansion.x + BOUND_OFFSET);
            maxCoordinate.x += (domain_expansion.y - BOUND_OFFSET);
            maxCoordinate.y += (domain_expansion.y - BOUND_OFFSET);
            if (dim == 3) {
                minCoordinate.z -= (domain_expansion.z + BOUND_OFFSET);
                maxCoordinate.z += (domain_expansion.z - BOUND_OFFSET);
            }
        }

        this->boundingBox[0] = minCoordinate;
        this->boundingBox[1] = maxCoordinate;

        generateGlobalBoundingBox(node_count);

        if (OPP_DBG)
            opp_printf("Local BoundingBox [computed]", 
                "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
                this->boundingBox[0].x, this->boundingBox[0].y, this->boundingBox[0].z, 
                this->boundingBox[1].x, this->boundingBox[1].y, this->boundingBox[1].z);
    }
    
    BoundingBox(const BoundingBox& other) : dim(other.dim) { OPP_RETURN_IF_INVALID_PROCESS;
        this->boundingBox[0] = other.boundingBox[0];
        this->boundingBox[1] = other.boundingBox[1];
        this->globalBoundingBox[0] = other.globalBoundingBox[0];
        this->globalBoundingBox[1] = other.globalBoundingBox[1];
        this->domain_expansion = other.domain_expansion;
    }

    //***********************************
    ~BoundingBox() { 
        // Nothing to implement
    }

    //***********************************
    inline const opp_point& getLocalMin() const {
        return this->boundingBox[0];
    }

    //***********************************
    inline const opp_point& getLocalMax() const {
        return this->boundingBox[1];
    }

    //***********************************
    inline const opp_point& getGlobalMin() const {
        return this->globalBoundingBox[0];
    }

    //***********************************
    inline const opp_point& getGlobalMax() const {
        return this->globalBoundingBox[1];
    }

    //***********************************
    inline bool isCoordinateInBoundingBox(const opp_point& point) {
        return !((point.x < boundingBox[0].x) || (point.x > boundingBox[1].x) ||
                (point.y < boundingBox[0].y) || (point.y > boundingBox[1].y) ||
                (point.z < boundingBox[0].z) || (point.z > boundingBox[1].z));
    }

    inline bool isCoordinateInGlobalBoundingBox(const opp_point& point) {
        return !((point.x < globalBoundingBox[0].x) || (point.x > globalBoundingBox[1].x) ||
                (point.y < globalBoundingBox[0].y) || (point.y > globalBoundingBox[1].y) ||
                (point.z < globalBoundingBox[0].z) || (point.z > globalBoundingBox[1].z));
    }

    //***********************************
    inline int getDim() const { 
        return dim; 
    }

private:
    //***********************************
    inline void generateGlobalBoundingBox(int count) {
#ifdef USE_MPI
        const double* localMin = reinterpret_cast<const double*>(&(this->boundingBox[0]));
        double* globalMin = reinterpret_cast<double*>(&(this->globalBoundingBox[0]));    
        MPI_Allreduce(localMin, globalMin, dim, MPI_DOUBLE, MPI_MIN, OPP_MPI_WORLD);

        const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
        double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
        MPI_Allreduce(localMax, globalMax, dim, MPI_DOUBLE, MPI_MAX, OPP_MPI_WORLD);

        // This is to avoid min and max corrdinates to not have MAX_REAL and MIN_REAL when current rank has no work
        if (count == 0) {
            this->boundingBox[0].x = this->globalBoundingBox[0].x; 
            this->boundingBox[0].y = this->globalBoundingBox[0].y; 
            this->boundingBox[0].z = this->globalBoundingBox[0].z; 

            this->boundingBox[1].x = this->globalBoundingBox[0].x;
            this->boundingBox[1].y = this->globalBoundingBox[0].y;
            this->boundingBox[1].z = this->globalBoundingBox[0].z;
        }

        if (OPP_rank == OPP_ROOT)
            opp_printf("Global BoundingBox", "Min[%2.6lE %2.6lE %2.6lE] Max[%2.6lE %2.6lE %2.6lE]", 
                this->globalBoundingBox[0].x, this->globalBoundingBox[0].y, this->globalBoundingBox[0].z, 
                this->globalBoundingBox[1].x, this->globalBoundingBox[1].y, this->globalBoundingBox[1].z);
#else
        this->globalBoundingBox[0] = this->boundingBox[0];
        this->globalBoundingBox[1] = this->boundingBox[1];
#endif
    }
};

//*******************************************************************************
class GlobalToLocalCellIndexMapper {

private:
    std::map<int,int> globalToLocalCellIndexMap;

public:
    //***********************************
    GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat, bool inc_halo = true) {
            
        globalToLocalCellIndexMap.clear();

        const opp_set c_set = global_cell_id_dat->set;
        const OPP_INT size = inc_halo ? 
                    (c_set->size + c_set->exec_size + c_set->nonexec_size) : (c_set->size);

        if (OPP_rank == 0)            
            opp_printf("GlobalToLocalCellIndexMapper", "Map Size %d", size);

        for (OPP_INT i = 0; i < size; i++) {
            const OPP_INT glbIndex = ((OPP_INT*)global_cell_id_dat->data)[i];
            globalToLocalCellIndexMap.insert(std::make_pair(glbIndex, i));
        }
        
        if (OPP_DBG && OPP_rank == 0)
            opp_printf("GlobalToLocalCellIndexMapper", "END");
    }

    //***********************************
    ~GlobalToLocalCellIndexMapper() {
        globalToLocalCellIndexMap.clear();

        if (OPP_DBG && OPP_rank == 0)
            opp_printf("~GlobalToLocalCellIndexMapper", "DONE");
    }

    //***********************************
    OPP_INT map(const OPP_INT globalIndex) const {
        if (globalIndex == MAX_CELL_INDEX)
            return MAX_CELL_INDEX;

        auto it = globalToLocalCellIndexMap.find(globalIndex);
        if (it != globalToLocalCellIndexMap.end())
            return it->second;
        
        if (OPP_DBG)
            opp_printf("GlobalToLocalCellIndexMapper", 
                "Error: Local cell index not found for global index %d", globalIndex);
        return MAX_INT;
    }
};

//*******************************************************************************
class CellMapper {

// private:
public:
    const std::shared_ptr<BoundingBox> boundingBox = nullptr;
    const double gridSpacing = 0.0;
    const double oneOverGridSpacing = 0.0;
    const opp_point& minGlbCoordinate;  
    const std::shared_ptr<Comm> comm = nullptr;
    const int dim;

    size_t globalGridDimsX = 0;
    size_t globalGridDimsY = 0;
    size_t globalGridDimsZ = 0;
    size_t globalGridDimsXY = 0;
    opp_uipoint localGridStart, localGridEnd;

    opp_uipoint minSavedCoordinate;
    opp_uipoint maxSavedCoordinate;

    size_t globalGridSize = 0;
    
    OPP_INT* structMeshToCellMapping = nullptr;         // This contain mapping to local cell indices
    OPP_INT* structMeshToRankMapping = nullptr;         // This contain mapping to residing mpi rank

    OPP_INT* structMeshToCellMapping_d = nullptr;       // This contain mapping to local cell indices
    OPP_INT* structMeshToRankMapping_d = nullptr;       // This contain mapping to residing mpi rank

    OPP_REAL* oneOverGridSpacing_d = nullptr;           // 1
    OPP_REAL* minGlbCoordinate_d = nullptr;             // DIM
    size_t* globalGridDims_d = nullptr;                 // DIM + 1 --> X, Y, Z, XY
    size_t* globalGridSize_d = nullptr;                 // 1

#ifdef USE_MPI        
    MPI_Win win_structMeshToCellMapping;
    MPI_Win win_structMeshToRankMapping;
#endif

public:
    // implementations in file opp_direct_hop_xxx.yyy
    CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, 
        const std::shared_ptr<Comm> comm = nullptr);
    ~CellMapper();

    void reduceInterNodeMappings(int callID);
    void createStructMeshMappingArrays();
    void convertToLocalMappings(const opp_dat global_cell_id_dat);
    void convertToLocalMappingsIncRank(const opp_dat global_cell_id_dat);
    void hostToDeviceTransferMappings();
    void generateStructuredMesh(opp_set set, const opp_dat c_gbl_id, 
            const std::function<void(const opp_point&, int&)>& all_cell_checker);
    void generateStructuredMeshFromFile(opp_set set, const opp_dat c_gbl_id);
    //***********************************
    inline void init_host(const double gridSpacing)
    {
        const opp_point& minGblCoordinate = boundingBox->getGlobalMin();
        const opp_point& maxGblCoordinate = boundingBox->getGlobalMax();

        size_t ax = 0, ay = 0, az = 0;

        { // removed this and added below due to decimal point issues
            // globalGridDimsX = std::ceil((maxGblCoordinate.x - minGblCoordinate.x) * oneOverGridSpacing);
            // globalGridDimsY = std::ceil((maxGblCoordinate.y - minGblCoordinate.y) * oneOverGridSpacing);
            // globalGridDimsZ = std::ceil((maxGblCoordinate.z - minGblCoordinate.z) * oneOverGridSpacing);   
        }
        {
            for (double z = minGblCoordinate.z; z < maxGblCoordinate.z; z += gridSpacing) az++;
            for (double y = minGblCoordinate.y; y < maxGblCoordinate.y; y += gridSpacing) ay++;
            for (double x = minGblCoordinate.x; x < maxGblCoordinate.x; x += gridSpacing) ax++; 
            globalGridDimsX = ax + 1;
            globalGridDimsY = ay + 1;
            globalGridDimsZ = az + 1; 
        }

        globalGridDimsXY = (globalGridDimsX * globalGridDimsY);
        globalGridSize = (globalGridDimsX * globalGridDimsY * globalGridDimsZ);

        if (OPP_rank == OPP_ROOT)
            opp_printf("CellMapper", "Global Grid Size - [%zu %zu %zu | %zu] gridSpacing [%2.10lE]", 
                globalGridDimsX, globalGridDimsY, globalGridDimsZ, globalGridDimsXY, gridSpacing); 
        
        const opp_point& minLocalCoordinate = boundingBox->getLocalMin();
        const opp_point& maxLocalCoordinate = boundingBox->getLocalMax();

        // Find the local ranks grid start indices
        ax = 0; ay = 0; az = 0;
        for (double z = minGblCoordinate.z; (z < minLocalCoordinate.z); z += gridSpacing) az++;
        for (double y = minGblCoordinate.y; (y < minLocalCoordinate.y); y += gridSpacing) ay++; 
        for (double x = minGblCoordinate.x; (x < minLocalCoordinate.x); x += gridSpacing) ax++; 
        localGridStart.x = ((ax == 0) ? 0 : (ax - 1));
        localGridStart.y = ((ay == 0) ? 0 : (ay - 1));
        localGridStart.z = ((az == 0) ? 0 : (az - 1));         

        // Find the local ranks grid end indices
        ax = 0; ay = 0; az = 0;
        for (double z = minGblCoordinate.z; (z <= maxLocalCoordinate.z); z += gridSpacing) az++; 
        for (double y = minGblCoordinate.y; (y <= maxLocalCoordinate.y); y += gridSpacing) ay++; 
        for (double x = minGblCoordinate.x; (x <= maxLocalCoordinate.x); x += gridSpacing) ax++; 
        localGridEnd.x = ((globalGridDimsX == ax) ? ax : (ax + 1));
        localGridEnd.y = ((globalGridDimsY == ay) ? ay : (ay + 1));
        localGridEnd.z = ((globalGridDimsZ == az) ? az : (az + 1));     

        const size_t local_grid_size = (size_t)(localGridEnd.x - localGridStart.x) * (size_t)(localGridEnd.y - localGridStart.y) * 
                                        (size_t)(localGridEnd.z - localGridStart.z);
        if (OPP_DBG)
            opp_printf("CellMapper", "Local Grid - Size [%zu] Min[%d %d %d] Max[%d %d %d]", 
                local_grid_size, localGridStart.x, localGridStart.y, localGridStart.z, 
                localGridEnd.x, localGridEnd.y, localGridEnd.z); 

        opp_minSavedDHGrid[0] = 0;
        opp_minSavedDHGrid[1] = 0;
        opp_minSavedDHGrid[2] = 0;

        opp_maxSavedDHGrid[0] = globalGridDimsX;
        opp_maxSavedDHGrid[1] = globalGridDimsY;
        opp_maxSavedDHGrid[2] = globalGridDimsZ;    
    }

    //***********************************
    inline opp_point getCentroidOfBox(const opp_point& coordinate) 
    { 
        opp_point centroid(MIN_REAL, MIN_REAL, MIN_REAL);
        const opp_point& maxCoordinate = boundingBox->getGlobalMax();

        switch (boundingBox->getDim()) {
            case 1:
                ASSIGN_CENTROID_TO_DIM(x); 
                break;
            case 2:
                ASSIGN_CENTROID_TO_DIM(x);
                ASSIGN_CENTROID_TO_DIM(y);
                break;
            case 3:
                ASSIGN_CENTROID_TO_DIM(x);
                ASSIGN_CENTROID_TO_DIM(y);
                ASSIGN_CENTROID_TO_DIM(z);
                break;
            default:
                std::cerr << "Error getCentroidOfBox: Dimension invalid " << 
                                boundingBox->getDim() << std::endl;
        }

        return centroid;
    }

    //***********************************
    inline size_t findStructuredCellIndex3D(const opp_point& position)  // Returns the global cell index
    { 
        // Round to the nearest integer to minimize rounding errors
        const size_t xIndex = static_cast<size_t>((position.x - minGlbCoordinate.x) * oneOverGridSpacing);
        const size_t yIndex = static_cast<size_t>((position.y - minGlbCoordinate.y) * oneOverGridSpacing);
        const size_t zIndex = static_cast<size_t>((position.z - minGlbCoordinate.z) * oneOverGridSpacing);

        bool isOutOfCuboid = (xIndex >= opp_maxSavedDHGrid[0] || xIndex < opp_minSavedDHGrid[0]) ||
                             (yIndex >= opp_maxSavedDHGrid[1] || yIndex < opp_minSavedDHGrid[1]) ||
                             (zIndex >= opp_maxSavedDHGrid[2] || zIndex < opp_minSavedDHGrid[2]);

        if (isOutOfCuboid) {
            return OPP_OUT_OF_SAVED_DOMAIN;
        }

        // Calculate the cell index mapping index
        const size_t index = xIndex + (yIndex * globalGridDimsX) + (zIndex * globalGridDimsXY);

        return (index >= globalGridSize) ? MAX_CELL_INDEX : index;
    }
    //***********************************
    inline size_t findStructuredCellIndex2D(const opp_point& position) // Returns the global cell index
    { 
        // Round to the nearest integer to minimize rounding errors
        const size_t xIndex = static_cast<size_t>((position.x - minGlbCoordinate.x) * oneOverGridSpacing);
        const size_t yIndex = static_cast<size_t>((position.y - minGlbCoordinate.y) * oneOverGridSpacing);

        // Calculate the cell index mapping index
        const size_t index = xIndex + (yIndex * globalGridDimsX);

        return (index >= globalGridSize) ? MAX_CELL_INDEX : index;
    }

    //***********************************
    inline int findClosestCellIndex(const size_t& structCellIdx) // Returns the global cell index
    { 
        if (OPP_DBG) {
            if (structCellIdx >= globalGridSize) {
                opp_printf("findClosestCellIndex", "Returning MAX - structCellIdx=%zu globalGridSize=%zu",
                    structCellIdx, globalGridSize);
                return MAX_CELL_INDEX;
            }
        }   
        return structMeshToCellMapping[structCellIdx];
    }

    //***********************************
    inline int findClosestCellRank(const size_t& structCellIdx) // Returns the rank of the cell
    { 
    #ifdef USE_MPI 
        if (OPP_DBG) {
            if (structCellIdx >= globalGridSize) {
                opp_printf("findClosestCellRank", "Returning MAX - structCellIdx=%zu globalGridSize=%zu",
                    structCellIdx, globalGridSize);
                return MAX_CELL_INDEX;
            }
        }
        return structMeshToRankMapping[structCellIdx];
    #else
        return OPP_rank;
    #endif
    }

    //***********************************
    inline void enrichStructuredMesh(const int index, const int cell_index, const int rank) 
    {
    #ifdef USE_MPI
        MPI_CHECK(MPI_Put(&cell_index, 1, MPI_INT, 0, index, 1, MPI_INT, win_structMeshToCellMapping));
        MPI_CHECK(MPI_Put(&rank, 1, MPI_INT, 0, index, 1, MPI_INT, win_structMeshToRankMapping));
    #else
        structMeshToCellMapping[index] = cell_index;
    #endif
    }

    //***********************************
    inline void printStructuredMesh(const std::string msg, int *array, size_t size, bool file = true, int line_break = 50) 
    {    
        if (!file) {
            opp_printf("structMeshToCellMapping", "%s - size=%zu", msg.c_str(), size);

            std::stringstream ss;
            for (size_t i = 0; i < size; i++) {
                ss << array[i] << ",";
                if ((i + 1) % line_break == 0) 
                    ss << "\n";
            }   
            printf("[RANK - %d]\n%s\n", OPP_rank, ss.str().c_str());
        }
        else {
            opp_write_array_to_file(array, size, msg);
        }
    }

    //***********************************
    inline void waitBarrier()
    {
    #ifdef USE_MPI
        MPI_CHECK(MPI_Barrier(OPP_MPI_WORLD));
        MPI_CHECK(MPI_Win_fence(0, win_structMeshToCellMapping)); 
        MPI_CHECK(MPI_Win_fence(0, win_structMeshToRankMapping)); 
    #endif
    }

    //***********************************
    inline void lockWindows() 
    {
    #ifdef USE_MPI
        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToCellMapping));
        MPI_CHECK(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win_structMeshToRankMapping));
    #endif
    }

    //***********************************
    inline void unlockWindows() 
    {
    #ifdef USE_MPI
        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToCellMapping));
        MPI_CHECK(MPI_Win_unlock(0, win_structMeshToRankMapping));
    #endif
    }

    inline void convertToLocalMappings_seq(const GlobalToLocalCellIndexMapper& mapper)
    {
        for (size_t i = 0; i < globalGridSize; i++) {
            if (structMeshToRankMapping[i] == OPP_rank) {   

                const int globalCID = (-1 * structMeshToCellMapping[i]);
                if ((globalCID != MAX_CELL_INDEX) || (globalCID != (-1 * MAX_CELL_INDEX))) {               
                    
                    const int localCID = mapper.map(globalCID);             
                    if (localCID != MAX_CELL_INDEX) {
                        structMeshToCellMapping[i] = localCID;
                    }
                    else {
                        opp_printf("CellMapper::convertToLocalMappings", 
                            "Error: cell mapping at %d is invalid [gcid:%d] rank_map %d rank %d", 
                            i, structMeshToCellMapping[i], structMeshToRankMapping[i], OPP_rank);
                    }
                }
            }
        }
    }
};

}; // end namespace opp