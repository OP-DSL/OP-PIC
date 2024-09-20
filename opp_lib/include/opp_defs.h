#pragma once

#include <limits>
#include <array>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <iostream>
#include <cstring>
#include <inttypes.h>
#include <stdio.h>
#include <stdarg.h>
#include <memory>
#include <cmath>
#include <sys/queue.h>

#ifdef USE_MPI
    namespace opp {
        class Comm;
        class GlobalParticleMover;
    };
#else
    #define Comm int
    #define GlobalParticleMover int
#endif

constexpr double MAX_REAL = std::numeric_limits<double>::max();
constexpr double MIN_REAL = std::numeric_limits<double>::min();

constexpr int MAX_INT = std::numeric_limits<int>::max();
constexpr int MIN_INT = std::numeric_limits<int>::min();

#define GET_VERT(D,K) ((K > maxCoordinate.D) ? maxCoordinate.D : K)

#define OPP_RUN_ON_ROOT(command) if (OPP_rank == OPP_ROOT command)

#define ASSIGN_CENTROID_TO_DIM(K)                                   \
    if (coordinate.K + this->gridSpacing <= maxCoordinate.K) {      \
        centroid.K = coordinate.K + this->gridSpacing * 0.5;         \
    }                                                               \
    else {                                                          \
        centroid.K = (coordinate.K + maxCoordinate.K) * 0.5;        \
    }                                                               \

#define CHECK(cmd) { int err = cmd; if (err != MPI_SUCCESS) opp_abort(std::to_string(err)); }

struct opp_point {
    opp_point(double _x, double _y, double _z) {
        x = _x; 
        y = _y;
        z = _z;
    };
    opp_point() { };

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct opp_ipoint {
    opp_ipoint(int _x, int _y, int _z) {
        x = _x; 
        y = _y;
        z = _z;
    };
    opp_ipoint() { };

    int x = 0;
    int y = 0;
    int z = 0;
};
