
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

#include <chrono>
#include <string>
#include <map>
#include <stack>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

// #define USE_MPI

#ifdef USE_MPI
#include <mpi.h>
extern MPI_Comm OPP_MPI_WORLD;
#endif

extern int OPP_rank;
extern int OPP_comm_size;

namespace opp {

    enum SetType {
        OPP_Mesh = 0,
        OPP_Particle,
    };

    struct ProfilerData {
        size_t count = 0;
        double value = 0.0;
    };

    class Profiler 
    {
    public:
        Profiler() {
            m_myRank = OPP_rank;
            m_worldSize = OPP_comm_size;
        }
        
        ~Profiler() {} 

        /**
        *  If it is unsure whether routine will be called on all MPI ranks, 
        *  register the profName prior to avoid deadlocks
        *  @param 
        */ 
        inline void reg(std::string profName) {
            // ProfilerData profData = m_elapsedTimeMap[profName];  

            m_elapsedTimeMap[profName] = ProfilerData();
            // m_elapsedTimeMap.insert(std::make_pair(profName, ProfilerData()));
            // Do nothing
        }

        /**
        *  
        *  @param 
        */        
        inline void start(std::string profName) { OPP_RETURN_IF_INVALID_PROCESS;
            m_startTimes[profName] = std::chrono::high_resolution_clock::now();
            m_currentProfileName.push(profName);
        }  

        /**
        *  
        *  @param 
        */ 
        inline double end(std::string profName) { OPP_RETURN_ZERO_IF_INVALID_PROCESS;
            if (profName == "")
                profName = m_currentProfileName.top();
#ifdef OPP_DEBUG
            if (m_currentProfileName.top() != profName)
                std::cerr << "opp::Profiler Error... Ending invalid profile ['" << profName << "' and '" << 
                    m_currentProfileName.top() << "']" << std::endl;
#endif
            m_currentProfileName.pop();
            double elapsedTime = this->getCurrentElapsedTime(profName, m_startTimes);
            ProfilerData& profData = m_elapsedTimeMap[profName];    
            profData.value += elapsedTime;
            profData.count += 1;
            return elapsedTime;
        } 

        /**
        *  
        *  @param 
        */ 
        inline std::string getCurrentProfileName() const {
            return m_currentProfileName.top();
        }

        /**
        *  
        *  @param 
        */ 
        inline void addTransferSize(std::string profName, SetType type, double transfer, size_t count = 0) {
            if (profName == "")
                profName = m_currentProfileName.top();
            ProfilerData* profData = nullptr;
            if (type == OPP_Mesh) {
                profData = &(m_meshTransferSizeMap[profName]);
            }
            else {
                profData = &(m_particleTransferSizeMap[profName]);
            }
            profData->value += transfer;
            profData->count += count;
        }

        /**
        *  
        *  @param 
        */ 
        inline void startMpiComm(std::string profName, SetType type) {
            if (profName == "")
                profName = m_currentProfileName.top();
            m_mpiCommStartTimes[profName + std::to_string(type)] = std::chrono::high_resolution_clock::now();
        }

        /**
        *  
        *  @param 
        */ 
        inline void endMpiComm(std::string profName, SetType type) {
            if (profName == "")
                profName = m_currentProfileName.top();
            ProfilerData* profData = nullptr;
            if (type == OPP_Mesh) {
                profData = &(m_meshMPITimeMap[profName]);
            }
            else {
                profData = &(m_particleMPITimeMap[profName]);
            }
            profData->value += getCurrentElapsedTime(profName + std::to_string(type), m_mpiCommStartTimes);
            profData->count += 1;
        }

        /**
        *  
        *  @param 
        */ 
        inline void printProfile(bool fromAllRanks = false) {
#ifndef USE_MPI
            fromAllRanks = true;
#endif
            m_myRank = OPP_rank;
            m_worldSize = OPP_comm_size;
            
            this->print(std::cout, fromAllRanks);
        }

        /**
        *  
        *  @param 
        */ 
        inline void printProfileToFile(std::string fileName, bool fromAllRanks = false) const {
            if (fromAllRanks)
                fileName += std::string("_") + std::to_string(m_myRank);
#ifndef USE_MPI
            fromAllRanks = true;
#endif            
            std::ofstream file(fileName);
            if (file.is_open()) {
                this->print(file, fromAllRanks);
                file.close();
            } else {
                std::cerr << "opp::Profiler Error... Failed to open file '" << fileName << "' for writing." << std::endl;
            }
        }

    private:
        inline double getCurrentElapsedTime(std::string profName, const std::map<std::string, std::chrono::high_resolution_clock::time_point>& start) const {
            std::chrono::duration<double> elapsedTime = (std::chrono::high_resolution_clock::now() - start.at(profName));
            return (double)elapsedTime.count();
        }

        inline double getElapsedTime(std::string profName) const {
            auto it = m_elapsedTimeMap.find(profName);
            if (it == m_elapsedTimeMap.end())
            {
                std::cerr << "opp::Profiler Error... profName: '" << profName << "' not found in rank " << m_myRank << std::endl;
                return 0.0;
            }

            return it->second.value;
        }

        inline void print(std::ostream& stream, bool fromAllRanks) const {  
            
            std::stringstream ss;
            if (fromAllRanks || m_myRank == 0) { 
                // ss << std::setprecision(3);
                ss << "Rank [" << (fromAllRanks ? std::to_string(m_myRank) : std::string("FINAL")) << "]"; 
                ss << "--------------------------------------------------------" << std::endl;

                ss << "'Profile'                      'Count'\t'TotalTime'\t'AverageTime'\t'MeshTransferSize'\t'MeshCommunicationTime'";
                ss << "\t'ParticleTransferSize'\t'ParticleTransferCount'\t'ParticleCommunicationTime'" << std::endl;
            }

            for (const auto& entry : m_elapsedTimeMap) {
                const std::string& profName = entry.first;

                // if (m_myRank == 0) printf("PROFILER ... Evaluating profile [%s]\n", profName.c_str());

                const size_t count = this->getMaxCallCount(profName, fromAllRanks); //entry.second.count;

                const double elapsedTime = this->getMaxElapsedTime(profName, fromAllRanks);
                const double averageTime = this->getAverageTime(profName);
                const double meshTransferSize = this->getMeshTransferSize(profName, fromAllRanks);
                const double meshCommunicationTime = this->getMeshCommunicationTime(profName, fromAllRanks);
                const double particleTransferSize = this->getParticleTransferSize(profName, fromAllRanks);
                const size_t particleTransferCount= this->getParticleTransferCount(profName, fromAllRanks);
                const double particleCommunicationTime = this->getParticleCommunicationTime(profName, fromAllRanks);

                if (fromAllRanks || m_myRank == 0) { 
                    ss << fromAllRanks << " " << m_myRank << " " << this->adjustString(profName, 30) << "\t" << this->adjustString(std::to_string(count), 6) << "\t"; 
                    ss << std::fixed << std::setprecision(10) << elapsedTime << "\t" << averageTime; 
                    std::stringstream so;
                    so << "\t" << meshTransferSize << "\t" << meshCommunicationTime;
                    so << "\t" << particleTransferSize << "\t" << particleTransferCount<< "\t" << particleCommunicationTime << std::endl;
                    ss << so.str();
                }
            }

            stream  << ss.str();
        }

        inline size_t getMaxCallCount(const std::string& profName, bool fromAllRanks) const {        
            return getMax(profName, m_elapsedTimeMap, fromAllRanks).count;
        }

        inline double getMaxElapsedTime(const std::string& profName, bool fromAllRanks) const {        
            return getMax(profName, m_elapsedTimeMap, fromAllRanks).value;
        }

        inline double getAverageTime(const std::string& profName) const {
            return (getSum(profName, m_elapsedTimeMap, false).value / m_worldSize);
        }

        inline double getMeshTransferSize(const std::string& profName, bool fromAllRanks) const {
            return getSum(profName, m_meshTransferSizeMap, fromAllRanks).value;
        }

        inline double getMeshCommunicationTime(const std::string& profName, bool fromAllRanks) const { 
            return getMax(profName, m_meshMPITimeMap, fromAllRanks).value;
        }

        inline double getParticleTransferSize(const std::string& profName, bool fromAllRanks) const {
            return getSum(profName, m_particleTransferSizeMap, fromAllRanks).value;
        }

        inline double getParticleTransferCount(const std::string& profName, bool fromAllRanks) const {
            return getSum(profName, m_particleTransferSizeMap, fromAllRanks).count;
        }

        inline double getParticleCommunicationTime(const std::string& profName, bool fromAllRanks) const {
            return getMax(profName, m_particleMPITimeMap, fromAllRanks).value;
        }

        inline ProfilerData getMax(const std::string& profName, const std::map<std::string, ProfilerData>& map, bool fromAllRanks) const {
            ProfilerData max;
            auto it = map.find(profName);
            if (it != map.end())
                max = it->second;
#ifdef USE_MPI
            if (!fromAllRanks)
            {
                ProfilerData localMax = max;  
                // Using all reduce since profile print is called after simulation - No perf issue
                MPI_Allreduce(&localMax.value, &max.value, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_WORLD); 
                MPI_Allreduce(&localMax.count, &max.count, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, OPP_MPI_WORLD);
            }
#endif  
            return max;
        }

        inline ProfilerData getSum(const std::string& profName, const std::map<std::string, ProfilerData>& map, bool fromAllRanks) const {
            ProfilerData sum;
            auto it = map.find(profName);
            if (it != map.end())
                sum = it->second;
#ifdef USE_MPI
            if (!fromAllRanks)
            {
                ProfilerData localSum = sum;  
                // Using all reduce since profile print is called after simulation - No perf issue
                MPI_Allreduce(&localSum.value, &sum.value, 1, MPI_DOUBLE, MPI_SUM, OPP_MPI_WORLD);
                MPI_Allreduce(&localSum.count, &sum.count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, OPP_MPI_WORLD);
            }
#endif  
            return sum;
        }

        inline std::string adjustString(const std::string& input, int maxLength) const {
            if (input.length() >= (size_t)maxLength) {
                return input.substr(0, maxLength);
            } else {
                std::string paddedString = input;
                paddedString.resize(maxLength, ' ');
                return paddedString;
            }
        }

        std::stack<std::string> m_currentProfileName;
        std::map<std::string, std::chrono::high_resolution_clock::time_point> m_startTimes;
        std::map<std::string, std::chrono::high_resolution_clock::time_point> m_mpiCommStartTimes;
        std::map<std::string, ProfilerData> m_elapsedTimeMap;
        std::map<std::string, ProfilerData> m_meshTransferSizeMap;
        std::map<std::string, ProfilerData> m_particleTransferSizeMap;
        std::map<std::string, ProfilerData> m_meshMPITimeMap;
        std::map<std::string, ProfilerData> m_particleMPITimeMap;

        int m_myRank = 0;
        int m_worldSize = 1;
    };

};

// int main(int argc, char **argv) {

// #ifdef USE_MPI
//     MPI_Init(&argc, &argv);    
// #endif

//     opp::Profiler timer;

//     timer.start("TIMER_1");
//     size_t counter = 0;
//     for (int i = 0; i < 1000000; i++)
//     {
//         if (i % 100000 == 0)
//         {
//             timer.start("TIMER_2");
//             for (int i = 0; i < 1000000; i++)    
//                 counter -= 1;
//             timer.end("TIMER_2", 10.1, 2.1, 0.01, 0.02);
//         }
//         counter += 1;
//     }
//     timer.end("TIMER_1", 10010.1, 1002.1, 1000.01, 1000.02);

//     timer.startTimer("TIMER_3");
//     counter = 0;
//     for (int i = 0; i < 1000000; i++)
//     {
//         if (i % 100000 == 0)
//         {
//             timer.startTimer("TIMER_4");
//             for (int i = 0; i < 1000000; i++)    
//                 counter -= 1;
//             timer.endTimer("TIMER_4");
//         }
//         counter += 1;
//     }
//     timer.endTimer("TIMER_3");

//     timer.printProfile(true);
//     timer.printProfileToFile("Rank", true);
//     timer.printProfileToFile("All", false);
// #ifdef USE_MPI
//     MPI_Finalize();
// #endif

//     return 0;
// }