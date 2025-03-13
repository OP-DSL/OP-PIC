/*==============================================================================*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 * Based on `fem-pic.cpp` by Lubos Brieda 
 * See https://www.particleincell.com/2015/fem-pic/ for more information
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#pragma once

#include <set>

/***************************************************************************************************
* maths.h
***************************************************************************************************/

/*compute inverse of a 3x3 matrix using the adjugate method*/
inline void inverse(double M[3][3], double V[3][3]) 
{

    double a=M[0][0];
    double b=M[0][1];
    double c=M[0][2];
    double d=M[1][0];
    double e=M[1][1];
    double f=M[1][2];
    double g=M[2][0];
    double h=M[2][1];
    double i=M[2][2];

    V[0][0]=(e*i-f*h);
    V[1][0]=-(d*i-f*g);
    V[2][0]=(d*h-e*g);
    V[0][1]=-(b*i-c*h);
    V[1][1]=(a*i-c*g);
    V[2][1]=-(a*h-b*g);
    V[0][2]=(b*f-c*e);
    V[1][2]=-(a*f-c*d);
    V[2][2]=(a*e-b*d);
    double det = a*V[0][0]+b*V[1][0]+c*V[2][0];

    double Vmax = 0;
    for (int m=0;  m<3; m++) 
    {
        for (int n=0;  n<3; n++) 
        {
            Vmax = fabs(V[m][n]) > Vmax ? fabs(V[m][n]) : Vmax;
        }
    }

    double idet=0;
    if (fabs(Vmax) / fabs(det) > 1e12) 
    {
        std::cerr<<"Matrix is not invertible, |det M| = " << fabs(det) << 
            "! setting to [0]."<<std::endl;
    }
    else 
        idet=1/det;

    /*1/det*/
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            V[i][j]*=idet;
}

/*computes determinant of a 3x3 matrix*/
inline double det3(double (*M)[3]) { 
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-
           M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+
           M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
}

inline double det4(double (*M)[4]) { 
    double M0[3][3];
    double M1[3][3];
    double M2[3][3];
    double M3[3][3];

    for (int i=0;i<3;i++) {
        M0[i][0]=M[i+1][1];
        M0[i][1]=M[i+1][2];
        M0[i][2]=M[i+1][3];

        M1[i][0]=M[i+1][0];
        M1[i][1]=M[i+1][2];
        M1[i][2]=M[i+1][3];

        M2[i][0]=M[i+1][0];
        M2[i][1]=M[i+1][1];
        M2[i][2]=M[i+1][3];

        M3[i][0]=M[i+1][0];
        M3[i][1]=M[i+1][1];
        M3[i][2]=M[i+1][2];
    }

    return M[0][0]*det3(M0) -
           M[0][1]*det3(M1) +
           M[0][2]*det3(M2) -
           M[0][3]*det3(M3);
}

/***************************************************************************************************
* particles.h
***************************************************************************************************/

/*particle*/
struct Particle {
    double pos[3];
    double vel[3];
    double lc[4];    /*particle's weights*/
    int cell_index;    /*last cell known to contain this particle*/
};

/***************************************************************************************************
* meshes.h
***************************************************************************************************/

//*********************************************
enum NodeType {
    NORMAL=0,
    OPEN=1,
    INLET=2,
    FIXED=3
};

//*********************************************
struct Node {
    Node(double x, double y, double z) {pos[0]=x;pos[1]=y;pos[2]=z;type=NORMAL;}
    double pos[3];    /*node position*/
    NodeType type;
    double volume = 0;    /*node volume*/
};

//*********************************************
struct Tetra {
    int con[4];
    double volume;
    Tetra (int n1, int n2, int n3, int n4) {con[0]=n1;con[1]=n2;con[2]=n3;con[3]=n4;}

    /*data structures to hold precomputed 3x3 determinants*/
    double alpha[4], beta[4], gamma[4], delta[4];

    /*cell connectivity*/
    int cell_con[4] = { -1, -1, -1, -1 };    /*index corresponds to the face opposite the i-th node*/
    bool initDone = false;
};

//*********************************************
struct Face {
    Face(int n1, int n2, int n3) {con[0]=n1, con[1]=n2, con[2]=n3;}
    int con[3];     // IDs of Nodes comprising the face
    double area;
    double u[3];
    double v[3];
    int cell_con = -1;
    double normal[3];
};

//*********************************************
struct Volume {
    std::vector <Node> nodes;
    std::vector <Tetra> elements;
    std::vector <Face> inlet_faces;
    double avg_edge_len;

    void summarize(std::ostream &out) {
        out << "MESH INFORMATION" << std::endl << "----------------" << std::endl;
        out << "  Number of nodes          = " << nodes.size() << std::endl;
        out << "  Number of elements       = " << elements.size() << std::endl;
        out << "  Number of inlet faces    = " << inlet_faces.size() << std::endl;
        out << std::endl;

        double min_pos[3] = {99999,99999,99999};
        double max_pos[3] = {0,0,0};
        for (auto node: nodes) {
            for (int d=0; d<3; d++) {
                min_pos[d] = std::min(node.pos[d], min_pos[d]);
                max_pos[d] = std::max(node.pos[d], max_pos[d]);
            }
        }

        out << "  Simulation region: "
            << max_pos[0] - min_pos[0] << "m x "
            << max_pos[1] - min_pos[1] << "m x "
            << max_pos[2] - min_pos[2] << "m" << std::endl;
        out << "  Average edge length: " << avg_edge_len <<"m"<< std::endl;

        out << std::endl;        
    }
};

//*********************************************
/*converts physical coordinate to logical
Returns true if particle matched to a tet
*/
inline bool XtoLtet(Particle &part, Volume &volume, bool search) {
    /*first try the current tetrahedron*/
    Tetra &tet = volume.elements[part.cell_index];

    bool inside = true;
    /*loop over vertices*/
    for (int i=0;i<4;i++) {
        part.lc[i] = (1.0/6.0)*(tet.alpha[i] - part.pos[0]*tet.beta[i] +
                      part.pos[1]*tet.gamma[i] - part.pos[2]*tet.delta[i])/tet.volume;
        if (part.lc[i]<0 || part.lc[i]>1.0) inside=false;
    }

    if (inside) return true;

    if (!search) return false;
    /*we are outside the last known tet, find most negative weight*/
    int min_i=0;
    double min_lc=part.lc[0];
    for (int i=1;i<4;i++)
        if (part.lc[i]<min_lc) {min_lc=part.lc[i];min_i=i;}

    /*is there a neighbor in this direction?*/
    if (tet.cell_con[min_i]>=0) {
        part.cell_index = tet.cell_con[min_i];
        return XtoLtet(part,volume,true);
    }

    return false;
}

//*********************************************
inline bool LoadVolumeMesh(const std::string file_name, Volume &volume) { 
    /*open file*/
    std::ifstream in(file_name);
    if (!in.is_open()) {std::cerr<<"Failed to open "<<file_name<<std::endl; return false;}

    printf("LoadVolumeMesh\n");

    /*read number of nodes and elements*/
    int n_nodes, n_elements;
    in>>n_nodes>>n_elements;

    /*read the nodes*/
    for (int n=0;n<n_nodes;n++) {
        int index;
        double x, y, z;

        in >> index >> x >> y >> z;
        if (index!=n+1) std::cerr<<"Inconsistent node numbering"<<std::endl;

        volume.nodes.emplace_back(x/1000.,y/1000.,z/1000.);

        if (n % 10000 == 0)
            printf("Read Nodes at %d\n", n);  
    }

    std::vector<double> edge_lengths;
    /*read elements, this will also contain edges and triangles*/
    for (int e=0;e<n_elements;e++) {
        int index, type;
        int n1, n2, n3, n4;

        in >> index >> type;

        if (type==304) { // Volume element
            in >> n1 >> n2 >> n3 >> n4;
            /*flipping nodes 2 & 3 to get positive volumes*/
            volume.elements.emplace_back(n1-1, n2-1, n3-1, n4-1);

        } 
        else {
            std::string s; getline(in,s);continue;
        }

        if (e % 10000 == 0)
            printf("Read Cells at %d\n", e);  
    }
    volume.avg_edge_len = 0;

    /*reset number of nodes and elements since we skipped bunch of lines and triangles*/
    n_nodes = volume.nodes.size();
    n_elements = volume.elements.size();

    std::map<int, std::vector<int>> node_con;

    /*compute element volumes*/
    for (int l=0; l<n_elements; l++) {

        Tetra &tet = volume.elements[l];

        {
            double M[4][4];

            /*set first column to 1*/
            for (int i=0;i<4;i++) M[i][0] = 1;

            /*loop over vertices*/
            for (int v=0;v<4;v++) {
                for (int dim=0;dim<3;dim++) {
                    M[v][dim+1] = volume.nodes[tet.con[v]].pos[dim];
                }

                node_con[tet.con[v]].emplace_back(l);
            }

            /*volume is (1/6)*det4(M)*/
            tet.volume = (1.0/6.0)*det4(M);

            /*flip ABCD to ADBC if negative volume*/
            if (tet.volume<0) {
                int t=tet.con[1];
                tet.con[1]=tet.con[3];
                tet.con[3]=t;
                tet.volume=-tet.volume;
            }
        }

        /*precompute 3x3 determinants for LC computation*/
        {
            double M[3][3];
            /*loop over vertices*/
            for (int v=0;v<4;v++) {
                int v2,v3,v4;

                switch (v) {
                    case 0: v2=1;v3=2;v4=3;break;
                    case 1: v2=3;v3=2;v4=0;break;
                    case 2: v2=3;v3=0;v4=1;break;
                    case 3: v2=1;v3=0;v4=2;break;
                }

                double *p2 = volume.nodes[tet.con[v2]].pos;
                double *p3 = volume.nodes[tet.con[v3]].pos;
                double *p4 = volume.nodes[tet.con[v4]].pos;

                /*alpha*/
                M[0][0] = p2[0];
                M[0][1] = p2[1];
                M[0][2] = p2[2];
                M[1][0] = p3[0];
                M[1][1] = p3[1];
                M[1][2] = p3[2];
                M[2][0] = p4[0];
                M[2][1] = p4[1];
                M[2][2] = p4[2];
                tet.alpha[v] = det3(M);

                /*beta*/
                M[0][0] =1;
                M[0][1] = p2[1];
                M[0][2] = p2[2];
                M[1][0] = 1;
                M[1][1] = p3[1];
                M[1][2] = p3[2];
                M[2][0] = 1;
                M[2][1] = p4[1];
                M[2][2] = p4[2];
                tet.beta[v] = det3(M);

                /*gamma*/
                M[0][0] =1;
                M[0][1] = p2[0];
                M[0][2] = p2[2];
                M[1][0] = 1;
                M[1][1] = p3[0];
                M[1][2] = p3[2];
                M[2][0] = 1;
                M[2][1] = p4[0];
                M[2][2] = p4[2];
                tet.gamma[v] = det3(M);

                /*delta*/
                M[0][0] =1;
                M[0][1] = p2[0];
                M[0][2] = p2[1];
                M[1][0] = 1;
                M[1][1] = p3[0];
                M[1][2] = p3[1];
                M[2][0] = 1;
                M[2][1] = p4[0];
                M[2][2] = p4[1];
                tet.delta[v] = det3(M);
            }
        }

        if (l % 100000 == 0)
            printf("Compute tetra at %d\n", l);  
    }

    for (int l=0;l<n_elements;l++) 
    {
        Tetra &tet = volume.elements[l];
        if (tet.initDone) continue;
        int v1,v2,v3;
        for (int v=0;v<4;v++) {
            /*skip if already set*/
            if (tet.cell_con[v]>=0) continue;

            switch(v) {
                case 0: v1=1;v2=2;v3=3;break;
                case 1: v1=2;v2=3;v3=0;break;
                case 2: v1=3;v2=0;v3=1;break;
                case 3: v1=0;v2=1;v3=2;break;
            }

            /*loop over the tets again looking for one with these three vertices*/
            /*(ejh - only look at elements that share a node with this element)*/
            for (int n=0; n<4; n++) {
                for (int m : node_con[tet.con[n]]) {
                    if (l == m) continue;
                    Tetra &other = volume.elements[m];
                    if (other.initDone) continue;

                    bool matches[4] = {false,false,false,false};
                    int count = 0;
                    for (int k=0;k<4;k++) {
                        if (other.con[k]==tet.con[v1] ||
                            other.con[k]==tet.con[v2] ||
                            other.con[k]==tet.con[v3]) {
                                count++;
                                matches[k]=true;
                        }
                    }

                    /*if three vertices match*/
                    if (count==3) {
                        tet.cell_con[v] = m;
                        if (tet.cell_con[0]>=0 && tet.cell_con[1]>=0 && tet.cell_con[2]>=0 && tet.cell_con[3]>=0)
                            tet.initDone = true;

                        /*set the cell connectivity for the index without a matching vertex to l*/
                        for (int k=0;k<4;k++)
                            if(!matches[k]) {
                                other.cell_con[k] = l;
                                if (other.cell_con[0]>=0 && other.cell_con[1]>=0 && other.cell_con[2]>=0 && other.cell_con[3]>=0)
                                    other.initDone = true;
                            }
                    }
                }
            }
        }

        if (l % 100000 == 0)
            printf("Compute tetra 2 at %d\n", l); 
    }

    /*also compute node volumes by scattering cell volumes,this can only be done after 3x3 dets are computed*/
    for (int i=0;i<n_elements;i++) {
        Particle dummy_part;
        Tetra &tet = volume.elements[i];
        dummy_part.cell_index = i;
        /*compute centroid position*/
        for (int dim=0;dim<3;dim++) {
            dummy_part.pos[dim]=0;
            for (int v=0;v<4;v++) dummy_part.pos[dim]+=0.25*volume.nodes[tet.con[v]].pos[dim];
        }

        bool found = XtoLtet(dummy_part,volume,false);
        if (!found) std::cerr<<"something is wrong"<<std::endl;

        for (int v=0;v<4;v++) {
            volume.nodes[tet.con[v]].volume += dummy_part.lc[v]*tet.volume;
        }

        /*mark nodes on open faces as open*/
        for (int v=0;v<4;v++) {
            if (tet.cell_con[v]<0)    /*no neighbor*/ {
                for (int i=0;i<4;i++) {
                    if (i!=v) volume.nodes[tet.con[i]].type=OPEN;
                }
            }
        }

        if (i % 100000 == 0)
            printf("Compute tetra 3 at %d\n", i); 
    }

    return true;
}

//*********************************************
/*loads nodes from a surface mesh file and sets them to the specified node type*/
inline bool LoadSurfaceMesh(const std::string file_name, Volume &volume, NodeType node_type, bool invert_normal) { 
    /*open file*/
    std::ifstream in(file_name);
    if (!in.is_open()) {std::cerr<<"Failed to open "<<file_name<<std::endl; return false;}

    printf("LoadSurfaceMesh\n"); 

    /*read number of nodes and elements*/
    int n_nodes, n_elements;
    in>>n_nodes>>n_elements;

    int nn = volume.nodes.size();

    /*read the nodes*/
    for (int n=0;n<n_nodes;n++) {
        int index;
        double x, y, z;

        in >> index >> x >> y >> z;

        if (index<1 || index>nn) {std::cerr<<"Incorrect node number "<<index<<std::endl;continue;}
        volume.nodes[index-1].type=node_type;

        if (n % 10000 == 0)
            printf("Read Surface Mesh at %d\n", n); 
    }

    if (node_type == INLET) {
        for (int e=0;e<n_elements;e++) {
            int index, type;
            int n1, n2, n3;

            in >> index >> type;

            if (type!=203) {std::string s; getline(in,s);continue;}

            in >> n1 >> n2 >> n3;

            /*flipping nodes 2 & 3 to get positive volumes*/
            volume.inlet_faces.emplace_back(n1-1, n2-1, n3-1);
            auto& inletFace = volume.inlet_faces.back();
            std::set<int> inletNodesSet = {n1 - 1, n2 - 1, n3 - 1};

            // Find the volume element that attaches to the inlet surface
            for (size_t v=0;v<volume.elements.size(); v++) {
                int matching_nodes = 0;
                for (int element_node: volume.elements[v].con) {
                    // if (element_node == n1-1 || element_node == n2-1 || element_node == n3-1) {
                    if (inletNodesSet.count(element_node)) {
                        matching_nodes += 1;
                        if (matching_nodes == 3) break;
                    }
                }
                if (matching_nodes == 3) {
                    if (inletFace.cell_con == -1) {
                         inletFace.cell_con = v;
                         break;
                    } else {
                        std::cerr<<"Inlet surface attached to more than one volume element"<<index<<std::endl;
                        exit(-1);
                    }
                }
            }
            if ( inletFace.cell_con == -1) {
                std::cerr<<"No volume element attached to inlet surface"<<index<<std::endl;
                exit(-1);
            }

            // Set the inlet velocity normal to the inlet surface
            for (int i=0; i<3; i++) {
                inletFace.u[i] = volume.nodes[n2-1].pos[i] - volume.nodes[n1-1].pos[i];
                inletFace.v[i] = volume.nodes[n3-1].pos[i] - volume.nodes[n1-1].pos[i];
            }

            double normal[3] = {
                inletFace.u[2] * inletFace.v[1] - inletFace.u[1] * inletFace.v[2],
                inletFace.u[0] * inletFace.v[2] - inletFace.u[2] * inletFace.v[0],
                inletFace.u[1] * inletFace.v[0] - inletFace.u[0] * inletFace.v[1]
            };

            // double normal[3];
            // normal[0] = inletFace.u[2]*inletFace.v[1] 
            //               - inletFace.u[1]*inletFace.v[2];
            // normal[1] = inletFace.u[0]*inletFace.v[2] 
            //               - inletFace.u[2]*inletFace.v[0];
            // normal[2] = inletFace.u[1]*inletFace.v[0] 
            //               - inletFace.u[0]*inletFace.v[1];
            if (invert_normal) {
                normal[0] = -normal[0];
                normal[1] = -normal[1];
                normal[2] = -normal[2];
            }

            double normal_len = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

            for (int i=0; i<3; i++) {
                inletFace.normal[i] = normal[i]/normal_len;
            }
            inletFace.area = normal_len / 2;

            if (e % 10000 == 0)
                printf("Compute inlet cells at %d\n", e); 
        }
    }

    return true;
}

//***************************************************************************************************