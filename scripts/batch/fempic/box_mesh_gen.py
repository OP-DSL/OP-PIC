import numpy as np
import time
# import pyvista as pv

# 102 : Edges
# 203 : Faces (Triangles)
# 304 : Tetrahedra

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Function '{func.__name__}' took {end_time - start_time} seconds to execute.")
        return result
    return wrapper

# --------------------------------------------------------
def get_node_index(x, y, z, nx, ny, nz):
    return (x) + (y) * (nx + 1) + (z) * (ny + 1) * (nx + 1)

# --------------------------------------------------------
@time_it
def create_nodes(nx, ny, nz, cube_length):

    # Calculate the total number of nodes = cube edges + centroids
    num_nodes = ((nx + 1) * (ny + 1) * (nz + 1)) + (nx * ny * nz)
    
    # Create node coordinates
    nodes = np.empty((num_nodes, 4), dtype=float)
    node_index = 0
    
    # Add corner vertices
    for z in range(nz + 1):
        for y in range(ny + 1):
            for x in range(nx + 1):
                nodes[node_index] = [(node_index+1), x * cube_length, y * cube_length, z * cube_length]
                node_index += 1
                
    # Calculate center coordinates for each cube and add them as nodes
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                center_x = x + 0.5
                center_y = y + 0.5
                center_z = z + 0.5
                center_node_index = (nx + 1) * (ny + 1) * (nz + 1) + z * nx * ny + y * nx + x
                nodes[center_node_index] = [(center_node_index+1), center_x * cube_length, center_y * cube_length, center_z * cube_length]

    return nodes

# --------------------------------------------------------
# Create tetrahedral connectivity centered around the center node
@time_it
def create_tetrahedra(nx, ny, nz):
    
    num_tetrahedra = nx * ny * nz * 12
    
    tetrahedra = np.empty((num_tetrahedra, 6), dtype=int)
    tetrahedra_index = 0
    
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):

                #define the 8 nodes bounding the cube
                p1 = get_node_index(x,   y,   z,    nx, ny, nz) + 1
                p2 = get_node_index(x,   y,   z+1,  nx, ny, nz) + 1
                p3 = get_node_index(x,   y+1, z,    nx, ny, nz) + 1
                p4 = get_node_index(x,   y+1, z+1,  nx, ny, nz) + 1
                p5 = get_node_index(x+1, y,   z,    nx, ny, nz) + 1
                p6 = get_node_index(x+1, y,   z+1,  nx, ny, nz) + 1
                p7 = get_node_index(x+1, y+1, z,    nx, ny, nz) + 1
                p8 = get_node_index(x+1, y+1, z+1,  nx, ny, nz) + 1
                
                center = ((nx + 1) * (ny + 1) * (nz + 1) + z * nx * ny + y * nx + x) +1

                # Define indices for the 12 tetrahedra #XY1 XY2, XY3, XY4, YZ1, YZ2, YZ3, YZ4,
                tetra = [[(tetrahedra_index+1), 304, p1, p3, p5, center],
                         [(tetrahedra_index+2), 304, p3, p5, p7, center],
                         [(tetrahedra_index+3), 304, p4, p2, p6, center],
                         [(tetrahedra_index+4), 304, p4, p6, p8, center],
                         [(tetrahedra_index+5), 304, p5, p6, p7, center],
                         [(tetrahedra_index+6), 304, p6, p7, p8, center],
                         [(tetrahedra_index+7), 304, p1, p3, p2, center],
                         [(tetrahedra_index+8), 304, p3, p4, p2, center],
                         [(tetrahedra_index+9), 304, p3, p4, p7, center],
                         [(tetrahedra_index+10), 304, p4, p7, p8, center],
                         [(tetrahedra_index+11), 304, p1, p5, p2, center],
                         [(tetrahedra_index+12), 304, p2, p5, p6, center]]

                for i in range(12):
                    tetrahedra[tetrahedra_index] = tetra[i]
                    tetrahedra_index += 1
    
    return tetrahedra

# --------------------------------------------------------
@time_it
def create_inlet_and_fixed_faces(nx, ny, nz, cube_length, nodes, tetrahedras):
    
    x = nx * cube_length
    y = ny * cube_length
    z = nz * cube_length

    boundary_nodes = []
    inlet_nodes = []
    inlet_faces = []
    
    # Seprate inlet and fixed boundary nodes
    for node in nodes:
        if (not(node[3] == 0 or node[3] == z)) and (node[1] == 0 or node[2] == 0 or node[1] == x or node[2] == y):
            boundary_nodes.append(node)
        
        if node[3] == 0:
            inlet_nodes.append(node)

    boundary_nodes = np.array(boundary_nodes, dtype=float)
    inlet_nodes = np.array(inlet_nodes, dtype=float)

    for inlet_node in inlet_nodes:

        matching_rows1 = tetrahedras[tetrahedras[:, 2] == inlet_node[0]]
        matching_rows2 = tetrahedras[tetrahedras[:, 3] == inlet_node[0]]
        matching_rows3 = tetrahedras[tetrahedras[:, 4] == inlet_node[0]]

        matching_rows = np.vstack((matching_rows1, matching_rows2, matching_rows3))

        for matching_row in matching_rows:
            # Mark as inlet face if all nodes are inlet nodes
            if np.any(inlet_nodes[:, 0] == matching_row[2]) and np.any(inlet_nodes[:, 0] == matching_row[3]) and np.any(inlet_nodes[:, 0] == matching_row[4]):
                inlet_faces.append(np.array([matching_row[2], matching_row[3], matching_row[4]]))
        
    inlet_faces = np.array(inlet_faces, dtype=np.int32)
    inlet_faces = np.unique(inlet_faces, axis=0)

    inlet_faces = np.hstack((np.full((inlet_faces.shape[0], 1), 4), np.full((inlet_faces.shape[0], 1), 203), inlet_faces))

    return boundary_nodes, inlet_nodes, inlet_faces

# --------------------------------------------------------
@time_it
def print_mesh(nodes, tetrahedras):
    file_path = "mesh.dat"
    with open(file_path, 'w') as file:
        
        file.write(f"{len(nodes)} {len(tetrahedras)}\n")

        for node in nodes:
            file.write(f"{int(node[0])} {'{:.14e}'.format(node[1])} {'{:.14e}'.format(node[2])} {'{:.14e}'.format(node[3])}\n")
        
        for tetrahedra in tetrahedras:
            file.write(f"{tetrahedra[0]} {tetrahedra[1]} {tetrahedra[2]} {tetrahedra[3]} {tetrahedra[4]} {tetrahedra[5]}\n")

# --------------------------------------------------------
@time_it
def print_inlet_faces(inlet_nodes, inlet_faces):
    file_path = "inlet.dat"
    with open(file_path, 'w') as file:
        
        file.write(f"{len(inlet_nodes)} {len(inlet_faces)}\n")

        for node in inlet_nodes:
            file.write(f"{int(node[0])} {'{:.14e}'.format(node[1])} {'{:.14e}'.format(node[2])} {'{:.14e}'.format(node[3])}\n")
        
        for face in inlet_faces:
            file.write(f"{face[0]} {face[1]} {face[2]} {face[3]} {face[4]}\n")   

# --------------------------------------------------------
@time_it
def print_fixed_nodes(boundary_nodes):
    file_path = "wall.dat"
    with open(file_path, 'w') as file:
        
        file.write(f"{len(boundary_nodes)} {0}\n")

        for node in boundary_nodes:
            file.write(f"{int(node[0])} {'{:.14e}'.format(node[1])} {'{:.14e}'.format(node[2])} {'{:.14e}'.format(node[3])}\n")
        
# --------------------------------------------------------
# 10 10 10 : 12000
# 10 20 10 : 24000
# 18 15 10 : 32400
# 20 20 10 : 48000
# 40 20 10 : 96000      1 - Node
# 40 40 10 : 192000     2 - Nodes
# 80 40 10 : 384000     4 - Nodes
# 80 80 10 : 768000     8 - Nodes
# 160 80 10 : 1536000   16 - Nodes
# 160 160 10 : 3072000  32 - Nodes
# 320 160 10 : 6144000  64 - Nodes
# 320 320 10 : 12288000  128 - Nodes
# 640 320 10 : 24576000  256 - Nodes
# 640 640 10 : 49152000  512 - Nodes

nx = 20                          # number of cubes in x direction
ny = 20*4                        # number of cubes in y direction
nz = 10                          # number of cubes in z direction
cube_length = 0.2                # in mm
tet_count = nx * ny * nz * 12    # there are 12 tetrahedrals per cube

nodes = create_nodes(nx, ny, nz, cube_length)
tetrahedras = create_tetrahedra(nx, ny, nz)
boundary_nodes, inlet_nodes, inlet_faces = create_inlet_and_fixed_faces(nx, ny, nz, cube_length, nodes, tetrahedras)

print_mesh(nodes, tetrahedras)
print_inlet_faces(inlet_nodes, inlet_faces)
print_fixed_nodes(boundary_nodes)

# --------------------------------------------------------

# celltypes = np.full(tet_count, pv.CellType.TETRA, dtype=np.uint8)
# pv_nodes = nodes[:, 1:]
# pv_tetrahedras = np.hstack((np.full((tetrahedras.shape[0], 1), 4), (tetrahedras[:, 2:] - 1)))

# mesh = pv.UnstructuredGrid(pv_tetrahedras, celltypes, pv_nodes)
# split_cells = mesh.explode(1)

# mesh.plot(show_edges=True, show_grid=True)
# split_cells.plot(show_edges=True, ssao=True, show_grid=True)

# print("nodes:" + str(len(nodes)) + "\n" + str(nodes))
# print("tetrahedras:" + str(len(tetrahedras)) + "\n" + str(tetrahedras))
# print("boundary_nodes:" + str(len(boundary_nodes)) + "\n" + str(boundary_nodes))
# print("inlet_nodes:" + str(len(inlet_nodes)) + "\n" + str(inlet_nodes))
# print("inlet_faces:" + str(len(inlet_faces)) + "\n" + str(inlet_faces))

# --------------------------------------------------------
