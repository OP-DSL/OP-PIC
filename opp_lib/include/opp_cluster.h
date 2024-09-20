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

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <unordered_map>
#include <algorithm> 

namespace opp {

struct Point3D {
    double x = 0.0f;
    double y = 0.0f;
    double z = 0.0f;

    bool operator==(const Point3D& other) const {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }
};

inline bool isEqual(const Point3D& p, const Point3D& other) {
    return (p.x == other.x) && (p.y == other.y) && (p.z == other.z);
}

// Function to calculate the euclidean distance
inline double euclideanDistance3D(const Point3D& p1, const Point3D& p2) {

    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate the euclidean distance
inline double euclideanDistancePlaneXY(const Point3D& p1, const Point3D& p2) {

    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Function to calculate the WCSS for a given clustering result
inline double calculateWCSS(const std::vector<Point3D>& data, 
    const std::vector<int>& clusterAssignments, const std::vector<Point3D>& centroids) {

    double wcss = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        int cluster = clusterAssignments[i];
        double distance = euclideanDistance3D(data[i], centroids[cluster]);
        wcss += distance * distance;
    }
    return wcss;
}


inline Point3D getTriangleCentroid3D(const Point3D& p1, const Point3D& p2, const Point3D& p3) {

    Point3D centroid;
    centroid.x = (p1.x + p2.x + p3.x) / 3.0;
    centroid.y = (p1.y + p2.y + p3.y) / 3.0;
    centroid.z = (p1.z + p2.z + p3.z) / 3.0;
    return centroid;
}

inline Point3D getTetraCentroid3D(const Point3D& p1, const Point3D& p2, const Point3D& p3, 
    const Point3D& p4) {

    Point3D centroid;
    centroid.x = (p1.x + p2.x + p3.x + p4.x) / 4.0;
    centroid.y = (p1.y + p2.y + p3.y + p4.y) / 4.0;
    centroid.z = (p1.z + p2.z + p3.z + p4.z) / 4.0;
    return centroid;
}

// Calculate the area of a triangle given its three vertices
inline double calculateTriangleArea(const Point3D& p1, const Point3D& p2, const Point3D& p3) {

    return 0.5 * ((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
}

// Check if a point is within a triangle using barycentric coordinates
inline bool isPointInTriangle(const Point3D& point, const Point3D& p1, const Point3D& p2, 
    const Point3D& p3) {

    double totalArea = calculateTriangleArea(p1, p2, p3);
    
    double alpha = calculateTriangleArea(point, p2, p3) / totalArea;
    double beta = calculateTriangleArea(p1, point, p3) / totalArea;
    double gamma = calculateTriangleArea(p1, p2, point) / totalArea;
    
    return (alpha >= 0 && beta >= 0 && gamma >= 0);
}

// Function to calculate the centroids of 3D points with cluster numbers
inline std::vector<Point3D> calculateTriangleCentroids3D(const std::vector<Point3D>& points, 
    const std::vector<int>& clusterNumbers) {

    std::unordered_map<int, std::vector<double>> centroids;
    std::unordered_map<int, int> clusterCounts;

    // Iterate over the points and update the centroids and cluster counts
    for (size_t i = 0; i < points.size(); ++i) {
        int cluster = clusterNumbers[i];
        const Point3D& point = points[i];

        centroids[cluster].resize(3, 0.0);  // Initialize centroid vector if it doesn't exist
        clusterCounts[cluster]++;  // Increment cluster count

        // Accumulate the coordinates for each dimension
        centroids[cluster][0] += point.x;
        centroids[cluster][1] += point.y;
        centroids[cluster][2] += point.z;
    }

    // Calculate the centroids by dividing the accumulated coordinates by the cluster counts
    std::vector<Point3D> result;
    for (const auto& entry : centroids) {
        int cluster = entry.first;
        const std::vector<double>& centroidSum = entry.second;
        int count = clusterCounts[cluster];

        Point3D point = { 
            (centroidSum[0] / count), 
            (centroidSum[1] / count), 
            (centroidSum[2] / count) 
        };

        result.push_back(point);
    }

    return result;
}

inline std::vector<int> BlockCluster(const std::vector<Point3D>& points, int numClusters) {

    std::vector<int> assignments; // Cluster assignments for each point

#ifdef USE_MPI
    int numPoint3Ds = (int)points.size();
    std::vector<int> clusterSizes(numClusters, 0);

    for (int r = 0; r < OPP_comm_size; r++) {

        int countForRank = opp_get_uniform_local_size(numPoint3Ds, r);

        for (int i = 0; i < countForRank; i++) {
            assignments.emplace_back(r);
            clusterSizes[r]++;
        }
    }
#endif

    // std::string logg = "";
    // for (auto a : clusterSizes) logg += std::to_string(a) + " ";
    // printf("BlockCluster SIZES %s", logg.c_str());

    return assignments;
}

inline std::vector<int> kMeansClustering3D(const std::vector<Point3D>& points, int numClusters) {
    int numPoint3Ds = points.size();

    // Randomly initialize cluster centroids
    std::vector<Point3D> centroids(numClusters);
    for (int i = 0; i < numClusters; ++i) {
        centroids[i] = points[std::rand() % numPoint3Ds];
    }

    std::vector<int> assignments(numPoint3Ds, 0); // Cluster assignments for each point
    std::vector<int> clusterSizes(numClusters, 0);
    std::vector<Point3D> newCentroids(numClusters, {0.0, 0.0, 0.0});

    bool converged = false;
    while (!converged) {

        std::fill(clusterSizes.begin(), clusterSizes.end(), 0);
        std::fill(newCentroids.begin(), newCentroids.end(), Point3D());

        // Assign points to clusters
        for (int i = 0; i < numPoint3Ds; ++i) {
            int closestCluster = 0;
            double minDistance = std::numeric_limits<double>::max();

            for (int j = 0; j < numClusters; ++j) {
                double dist = euclideanDistance3D(points[i], centroids[j]);
                if (dist < minDistance) {
                    minDistance = dist;
                    closestCluster = j;
                }
            }

            assignments[i] = closestCluster;
            clusterSizes[closestCluster]++;
            newCentroids[closestCluster].x += points[i].x;
            newCentroids[closestCluster].y += points[i].y;
            newCentroids[closestCluster].z += points[i].z;
        }

        // Update cluster centroids
        for (int i = 0; i < numClusters; ++i) {
            if (clusterSizes[i] > 0) {
                newCentroids[i].x /= clusterSizes[i];
                newCentroids[i].y /= clusterSizes[i];
                newCentroids[i].z /= clusterSizes[i];
            }
        }

        // Check convergence
        converged = (centroids == newCentroids);
        centroids = newCentroids;
    }

    std::string logg = "";
    for (auto a : clusterSizes) logg += std::to_string(a) + " ";
    printf("SIZES %s", logg.c_str());

    return assignments;
}

// Function to calculate the volume of a tetrahedron
inline double calculateTetraVolume(const Point3D& p1, const Point3D& p2,
                                  const Point3D& p3, const Point3D& p4) {
    // Calculate the vectors between the vertices
    double v1x = p2.x - p1.x;
    double v1y = p2.y - p1.y;
    double v1z = p2.z - p1.z;
    
    double v2x = p3.x - p1.x;
    double v2y = p3.y - p1.y;
    double v2z = p3.z - p1.z;
    
    double v3x = p4.x - p1.x;
    double v3y = p4.y - p1.y;
    double v3z = p4.z - p1.z;
    
    // Calculate the volume using the scalar triple product
    double volume = std::abs(v1x * (v2y * v3z - v3y * v2z)
                             - v1y * (v2x * v3z - v3x * v2z)
                             + v1z * (v2x * v3y - v3x * v2y)) / 6.0;
    
    return volume;
}

};