
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
    //***********************************
    BoundingBox(int dim, opp_point minCoordinate, opp_point maxCoordinate) : dim(dim) {

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
    BoundingBox(const opp_dat node_pos_dat, int dim, const OPP_REAL expansion[3]) : dim(dim) {
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
            minCoordinate.x -= (expansion[0] + BOUND_OFFSET);
            minCoordinate.y -= (expansion[0] + BOUND_OFFSET);
            maxCoordinate.x += (expansion[1] - BOUND_OFFSET);
            maxCoordinate.y += (expansion[1] - BOUND_OFFSET);
            if (dim == 3) {
                minCoordinate.z -= (expansion[2] + BOUND_OFFSET);
                maxCoordinate.z += (expansion[2] - BOUND_OFFSET);
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
    
    BoundingBox(const BoundingBox& other) : dim(other.dim) {
        this->boundingBox[0] = other.boundingBox[0];
        this->boundingBox[1] = other.boundingBox[1];
        this->globalBoundingBox[0] = other.globalBoundingBox[0];
        this->globalBoundingBox[1] = other.globalBoundingBox[1];
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
        MPI_Allreduce(localMin, globalMin, dim, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        const double* localMax = reinterpret_cast<const double*>(&(this->boundingBox[1]));
        double* globalMax = reinterpret_cast<double*>(&(this->globalBoundingBox[1]));    
        MPI_Allreduce(localMax, globalMax, dim, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

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
    GlobalToLocalCellIndexMapper(const opp_dat global_cell_id_dat) {
        if (OPP_rank == 0)            
            opp_printf("GlobalToLocalCellIndexMapper", "START");
            
        globalToLocalCellIndexMap.clear();

        const opp_set cells_set = global_cell_id_dat->set;
        const OPP_INT size_inc_halo = (cells_set->size + cells_set->exec_size + cells_set->nonexec_size);

        for (OPP_INT i = 0; i < size_inc_halo; i++) {
            const OPP_INT glbIndex = ((OPP_INT*)global_cell_id_dat->data)[i];
            globalToLocalCellIndexMap.insert(std::make_pair(glbIndex, i));
        }
        
        if (OPP_rank == 0)
            opp_printf("GlobalToLocalCellIndexMapper", "END");
    }

    //***********************************
    ~GlobalToLocalCellIndexMapper() {
        globalToLocalCellIndexMap.clear();

        if (OPP_rank == 0)
            opp_printf("~GlobalToLocalCellIndexMapper", "DONE");
    }

    //***********************************
    OPP_INT map(const OPP_INT globalIndex) {
        if (globalIndex == MAX_CELL_INDEX)
            return MAX_CELL_INDEX;

        auto it = globalToLocalCellIndexMap.find(globalIndex);
        if (it != globalToLocalCellIndexMap.end())
            return it->second;
        
        opp_printf("GlobalToLocalCellIndexMapper", 
            "Error: Local cell index not found for global index %d", globalIndex);
        return MAX_INT;
    }
};

//*******************************************************************************
class CellMapper {

public:
    // implementations in file opp_direct_hop_core.cpp
    CellMapper(const std::shared_ptr<BoundingBox> boundingBox, const double gridSpacing, 
        const std::shared_ptr<Comm> comm = nullptr);
    ~CellMapper();

    opp_point getCentroidOfBox(const opp_point& coordinate);
    void reduceInterNodeMappings(int callID);
    void convertToLocalMappings(const opp_dat global_cell_id_dat);
    void createStructMeshMappingArrays();

    //***********************************
    inline size_t findStructuredCellIndex3D(const opp_point& position)  // Returns the global cell index
    { 
        // Round to the nearest integer to minimize rounding errors
        const size_t xIndex = static_cast<size_t>((position.x - minGlbCoordinate.x) * oneOverGridSpacing);
        const size_t yIndex = static_cast<size_t>((position.y - minGlbCoordinate.y) * oneOverGridSpacing);
        const size_t zIndex = static_cast<size_t>((position.z - minGlbCoordinate.z) * oneOverGridSpacing);

        // Calculate the cell index mapping index
        const size_t index = xIndex + (yIndex * globalGridDimsX) + (zIndex * globalGridDimsXY);

        return (index >= globalGridSize || index < 0) ? MAX_CELL_INDEX : index;
    }
    //***********************************
    inline size_t findStructuredCellIndex2D(const opp_point& position) // Returns the global cell index
    { 
        // Round to the nearest integer to minimize rounding errors
        const int xIndex = static_cast<int>((position.x - minGlbCoordinate.x) * oneOverGridSpacing);
        const int yIndex = static_cast<int>((position.y - minGlbCoordinate.y) * oneOverGridSpacing);

        // Calculate the cell index mapping index
        const size_t index = ((size_t)(xIndex) + ((size_t)yIndex * globalGridDimsX));

        return (index >= globalGridSize || index < 0) ? MAX_CELL_INDEX : index;
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
            printf("%s\n", ss.str().c_str());
        }
        else {
            opp_write_array_to_file(array, size, msg);
        }
    }

    //***********************************
    inline void waitBarrier()
    {
    #ifdef USE_MPI
        MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));
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

// private:
    const std::shared_ptr<BoundingBox> boundingBox = nullptr;
    const double gridSpacing = 0.0;
    const double oneOverGridSpacing = 0.0;
    const opp_point& minGlbCoordinate;  
    const std::shared_ptr<Comm> comm = nullptr;

    size_t globalGridDimsX = 0;
    size_t globalGridDimsY = 0;
    size_t globalGridDimsZ = 0;
    size_t globalGridDimsXY = 0;
    opp_uipoint localGridStart, localGridEnd;

    size_t globalGridSize = 0;
    
    int* structMeshToCellMapping = nullptr;         // This contain mapping to local cell indices
    int* structMeshToRankMapping = nullptr;         // This contain mapping to residing mpi rank

#ifdef USE_MPI        
    MPI_Win win_structMeshToCellMapping;
    MPI_Win win_structMeshToRankMapping;
#endif
};

}; // end namespace opp