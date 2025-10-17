import sys

# This code assumes that the filtration is given in an ASCII file with the following format, where each line represents a simplex sigma and is of the form:
# f(sigma) dim(sigma) v_0 ... v_{dim(sigma)}
# where:
# f(sigma) is the function value of sigma (its "time of appearance" in the filtration),
# dim(sigma) is the dimension of sigma,
# v_0 ... v_{dim(sigma)} are the IDs (integers) of the vertices of sigma.
# For instance, 0.125 2 2 6 4 denotes a 2-simplex (triangle) that appears in the filtration at time 0.125 and whose vertices have IDs 2, 6 and 4. Warning: the function values provided in the file must be compatible with the underlying simplicial complex (function values of simplices must be at least as large as the ones of their faces). Nevertheless, the vertex IDs are arbitrary integers and may not start at 0 nor be continuous.

class Simplex:
    def __init__(self, tokens: list[str], index: int):
        """
        Read a simplex from tokens starting at index.
        Returns the next index to read from.
        """
        self.val = float(tokens[index])
        self.dim = int(tokens[index + 1])
        self.vert = set()
        
        for i in range(self.dim + 1):
            self.vert.add(int(tokens[index + 2 + i]))
        
        # Return next index
        self.next_index = index + 2 + self.dim + 1
    
    def __str__(self):
        return f"{{val={self.val}; dim={self.dim}; {sorted(self.vert)}}}\n"
    
    def __repr__(self):
        return self.__str__()

FilteredSimplicialComplex = list[Simplex]

def read_filtration(filename: str) -> FilteredSimplicialComplex:
    """Read filtration from file and return list of Simplex objects."""
    filtration = []
    
    with open(filename, 'r') as file:
        # Read all tokens (whitespace-separated values)
        content = file.read()
        tokens = content.split()
    
    # Parse tokens into Simplex objects
    index = 0
    while index < len(tokens):
        simplex = Simplex(tokens, index)
        filtration.append(simplex)
        index = simplex.next_index
    
    return filtration


# def main():
#     if len(sys.argv) != 2:
#         print("Syntax: python read_filtration.py <filename>")
#         sys.exit(0)
    
#     try:
#         filtration = read_filtration(sys.argv[1])
#         print(filtration)
#     except FileNotFoundError:
#         print(f"Error: File '{sys.argv[1]}' not found")
#         sys.exit(1)


# if __name__ == "__main__":
#     main()

def sort_simplices(filtration: FilteredSimplicialComplex):
    """Sort simplices in filtration by (val, lexicographic order of vertices)."""
    def simplex_key(simplex):
        return (simplex.val, sorted(simplex.vert))
    
    filtration.sort(key=simplex_key)

# TODO: better to insert a simplex that is going to be killed directly after
# Take this into account when sorting the indices

# Test sorting
# def main():
#     if len(sys.argv) != 2:
#         print("Syntax: python read_filtration.py <filename>")
#         sys.exit(0)
    
#     try:
#         filtration = read_filtration(sys.argv[1])
#         print("Before sorting:")
#         print(filtration)
        
#         sort_simplices(filtration)
        
#         print("After sorting:")
#         print(filtration)
#     except FileNotFoundError:
#         print(f"Error: File '{sys.argv[1]}' not found")
#         sys.exit(1)

# if __name__ == "__main__":
#     main()

SparseMatrix = list[set[int]]

# 1. Build boundary matrix
def build_boundary_matrix_sparse(filtration: FilteredSimplicialComplex) -> SparseMatrix:
    """Build the boundary matrix in sparse format from the filtration.
    Returns a list of columns, where each column is represented as a set of row indices with non-zero entries."""
    boundary_matrix = []
    
    simplex_lookup = {} # Store previous simplices for quick lookup, we suppose filtration is sorted

    for i, simplex in enumerate(filtration):
        simplex_lookup[(simplex.dim, str(sorted(simplex.vert)))] = i # need hashable key so no list or set
        # sorted ensures {1,2} and {2,1} are treated the same
        if simplex.dim == 0:
            # 0-simplices have no boundary
            boundary_matrix.append(set())
        else:
            # Compute boundary as set of row indices
            boundary = set()
            for v in simplex.vert:
                # Create face by removing vertex v
                face = simplex.vert - {v}
                # Find index of face in filtration
                # Before :
                # for i, s in enumerate(filtration):
                #     if s.dim == simplex.dim - 1 and s.vert == face:
                #         boundary.add(i)
                #         break
                # After:
                face_index = simplex_lookup.get((simplex.dim - 1, str(sorted(face))))
                if face_index is not None:
                    boundary.add(face_index)
            boundary_matrix.append(boundary)
    
    return boundary_matrix

DenseMatrix = list[list[int]]

# Utils: convert sparse matrix to dense matrix
def sparse_to_dense(matrix: SparseMatrix, n_rows: int) -> DenseMatrix:
    """Convert sparse matrix (list of sets) to dense matrix (list of lists i.e. list of rows)."""
    dense_matrix = []

    for i in range(n_rows):
        dense_matrix.append([0] * len(matrix))
        for j in range(len(matrix)):
            if i in matrix[j]:
                dense_matrix[i][j] = 1
    
    return dense_matrix

# Utils: log dense matrix
def log_dense_matrix(matrix: DenseMatrix):
    """Log the dense matrix to the console."""
    for row in matrix:
        print(" ".join(str(x) for x in row))

# Test building boundary matrix
# def main():
#     if len(sys.argv) != 2:
#         print("Syntax: python read_filtration.py <filename>")
#         sys.exit(0)
    
#     try:
#         filtration = read_filtration(sys.argv[1])
#         sort_simplices(filtration)
#         print(filtration)

#         boundary_matrix = build_boundary_matrix_sparse(filtration)
        
#         print("Boundary matrix (sparse format):")
#         for i, col in enumerate(boundary_matrix):
#             print(f"Column {i}: {sorted(col)}")
#     except FileNotFoundError:
#         print(f"Error: File '{sys.argv[1]}' not found")
#         sys.exit(1)

# if __name__ == "__main__":
#     main()

# 2. Reduction algorithm of the boundary matrix
def reduce_matrix_sparse(matrix: SparseMatrix) -> tuple[SparseMatrix, list[tuple[int, int]]]:
    """Reduce the boundary matrix in sparse format using column operations.
    Returns the reduced boundary matrix and a list of low indices for each column."""
    m = len(matrix)
    low = [] # List of couples (column index, low index) for each column that is non null
    
    for j in range(m): # Complexity : linear (number of entries), composed : n^3
        current_low = -1

        while True: # Complexity : if found each time, linear (number of entries), composed : n^2
            if not matrix[j]:
                break  # Column is zero, low does not contain the column, and we don't add this null column to it
            
            # Else column is non zero, find the lowest index Find lowest index in column j
            current_low = max(matrix[j])
            
            # while i < j st low(M_-,i) = low(M_-,j) and non zero
            found = False
            for column_index, low_index in low: # Complexity : linear (number of entries), composed : one execution of symmetric difference, so 2*number of entries
                if low_index == current_low:
                    # Add column i to column j (symmetric difference)
                    matrix[j] ^= matrix[column_index] # Complexity : linear (number of entries), executed once per loop at most (break after)
                    found = True
                    break

            if not found:
                low.append((j, current_low))
                break  # No more reductions possible
    
    return matrix, low

# Complexity analysis: O(n^3) in worst case
# FIXME: 
# - In the case of sparse matrices, consider symmetric difference is constant
# - For the i loop, since we count number of operations on columns happens only once per loop so still constant for this loop
# - Note that storing only low indices where the columns are not null helps reduce the number of iterations in practice and must be taken into account for sparse matrix
# - While True : remonte le pivot autant de fois qu'il y a de low, c'est constant
# - La outer boucle donne linéaire

# Test reduction algorithm
# def main():
#     if len(sys.argv) != 2:
#         print("Syntax: python read_filtration.py <filename>")
#         sys.exit(0)
    
#     try:
#         filtration = read_filtration(sys.argv[1])
#         sort_simplices(filtration)
#         print(filtration)

#         boundary_matrix = build_boundary_matrix_sparse(filtration)

#         dense_matrix = sparse_to_dense(boundary_matrix, len(filtration))
#         log_dense_matrix(dense_matrix)
        
#         print("Boundary matrix (sparse format):")
#         for i, col in enumerate(boundary_matrix):
#             print(f"Column {i}: {sorted(col)}")
        
#         reduced_matrix, low = reduce_matrix_sparse(boundary_matrix)

#         dense_reduced_matrix = sparse_to_dense(reduced_matrix, len(filtration))
#         log_dense_matrix(dense_reduced_matrix)
        
#         print("Reduced boundary matrix (sparse format):")
#         for i, col in enumerate(reduced_matrix):
#             print(f"Column {i}: {sorted(col)}")
        
#         print("Low indices:", low)
#     except FileNotFoundError:
#         print(f"Error: File '{sys.argv[1]}' not found")
#         sys.exit(1)

# if __name__ == "__main__":
#     main()

# 4. Output the barcode : 
# - 1 line per interval, containing 3 numbers: 
# - dimension of the feature, left endpoint (filtration value associated with the simplex that created the homological feature) and 
# - right endpoint (filtration value associated with the simple that killed the feature) separated by white spaces
def compute_barcode(matrix: SparseMatrix, low: list[tuple[int, int]], filtration: FilteredSimplicialComplex) -> list[tuple[int, int, int]]:
    """Compute the barcode from the reduced boundary matrix and low indices.
    Returns a list of intervals (dim, birth, death). The filtration is needed to get the dimension of the original simplices.
    It is supposed that the filtration is sorted."""
    barcode = []
    m = len(matrix)
    low_dict = {col: low_idx for col, low_idx in low}
    
    # Invariant in low : col > low_idx
    # The low index gives birth, the column kills
    genitors_killed = set(low_idx for _, low_idx in low) # killers are the keys in low_dict

    for j in range(m):
        simplex = filtration[j]
        if j in low_dict:
            # Column is non-zero, feature dies
            low_index = low_dict[j]
            birth_simplex = filtration[low_index]
            barcode.append((birth_simplex.dim, birth_simplex.val, simplex.val))
        elif j not in genitors_killed:
            # Column is zero, feature persists to infinity
            barcode.append((simplex.dim, simplex.val, float('inf')))
    
    return barcode

def write_barcode(barcode: list[tuple[int, int, int]], filename: str):
    """Write the barcode to a file."""
    with open(filename, 'w') as file:
        for dim, birth, death in barcode:
            file.write(f"{dim} {birth} {death}\n")

from matplotlib import pyplot as plt

def plot_barcode(barcode: list[tuple[int, int, int]]):
    """Plot the barcode using matplotlib.
    Death at infinite is plotted with a red arrow."""
    plt.figure()
    y = 0
    for dim, birth, death in barcode:
        if death != float('inf'):
            plt.hlines(y, birth, death, colors='b', lw=2)
        else:
            plt.hlines(y, birth, birth + 1, colors='b', lw=2)
            plt.arrow(birth + 1, y, 0.5, 0, head_width=1, head_length=0.1, fc='r', ec='r')
        y += 1
    plt.xlabel('Filtration Value')
    plt.ylabel('Homological Features')
    plt.title('Barcode')
    plt.show()

# def main():
#     if len(sys.argv) != 3:
#         print("Syntax: python main.py <input_filename> <output_filename>")
#         sys.exit(0)
    
#     try:
#         filtration = read_filtration(sys.argv[1])
#         sort_simplices(filtration)

#         boundary_matrix = build_boundary_matrix_sparse(filtration)
        
#         reduced_matrix, low = reduce_matrix_sparse(boundary_matrix)

#         barcode = compute_barcode(reduced_matrix, low, filtration)

#         write_barcode(barcode, sys.argv[2])
        
#         print(f"Barcode written to '{sys.argv[2]}'")
#     except FileNotFoundError:
#         print(f"Error: File '{sys.argv[1]}' not found")
#         sys.exit(1)

# if __name__ == "__main__":
#     main()

# python .\main.py .\filtrations\triangle .\output\triangle

# 5. d-sphere, moebius band, torus, klein bottle and projective plane

from itertools import combinations

def build_sphere(d: int) -> list[tuple[float, int, list[int]]]:
    """Build a (d-1)-sphere as the boundary of a d-simplex."""
    complex = []
    vertices = list(range(d + 1))
    
    # Generate all faces of dimension 0 to d-1, dim d is for the disk
    for dim in range(d):
        for simplex_vertices in combinations(vertices, dim + 1):
            complex.append((float(dim), dim, list(simplex_vertices)))
    
    return complex

def build_torus():
    """Build a torus"""
    complex = []

    # Form 3 triangles
    complex.append((0.0, 0, [0]))
    complex.append((1.0, 0, [1]))
    complex.append((2.0, 0, [2]))
    complex.append((3.0, 1, [0, 1]))
    complex.append((4.0, 1, [1, 2]))
    complex.append((5.0, 1, [2, 0]))
    complex.append((6.0, 0, [3]))
    complex.append((7.0, 0, [4]))
    complex.append((8.0, 0, [5]))
    complex.append((9.0, 1, [3, 4]))
    complex.append((10.0, 1, [4, 5]))
    complex.append((11.0, 1, [5, 3]))
    complex.append((12.0, 0, [6]))
    complex.append((13.0, 0, [7]))
    complex.append((14.0, 0, [8]))
    complex.append((15.0, 1, [6, 7]))
    complex.append((16.0, 1, [7, 8]))
    complex.append((17.0, 1, [8, 6]))

    # FIXME: I may have added twists under there
    # Link the triangles with twisted edges
    # Between triangle 1 and 2
    complex.append((18.0, 1, [0, 3]))
    complex.append((19.0, 1, [1, 4]))
    complex.append((20.0, 1, [2, 5]))
    # Add the faces
    complex.append((21.0, 2, [0, 1, 4])) # two edges from triangle 1 and one from triangle 2
    complex.append((22.0, 2, [1, 2, 5]))
    complex.append((23.0, 2, [2, 0, 3]))
    complex.append((24.0, 2, [0, 4, 3])) # two edges from triangle 2 and one from triangle 1
    complex.append((25.0, 2, [1, 5, 4]))
    complex.append((26.0, 2, [2, 3, 5]))

    # Between triangle 2 and 3
    complex.append((27.0, 1, [3, 6]))
    complex.append((28.0, 1, [4, 7]))
    complex.append((29.0, 1, [5, 8]))
    # Add the faces
    complex.append((30.0, 2, [3, 4, 7])) # two edges from triangle 2 and one from triangle 3
    complex.append((31.0, 2, [4, 5, 8]))
    complex.append((32.0, 2, [5, 3, 6]))
    complex.append((33.0, 2, [3, 7, 6])) # two edges from triangle 3 and one from triangle 2
    complex.append((34.0, 2, [4, 8, 7]))
    complex.append((35.0, 2, [5, 6, 8]))

    # Between triangle 3 and 1
    complex.append((36.0, 1, [6, 0]))
    complex.append((37.0, 1, [7, 1]))
    complex.append((38.0, 1, [8, 2]))
    # Add the faces
    complex.append((39.0, 2, [6, 7, 1])) # two edges from triangle 3 and one from triangle 1
    complex.append((40.0, 2, [7, 8, 2]))
    complex.append((41.0, 2, [8, 6, 0]))
    complex.append((42.0, 2, [6, 1, 0])) # two edges from triangle 1 and one from triangle 3
    complex.append((43.0, 2, [7, 2, 1]))
    complex.append((44.0, 2, [8, 0, 2]))

    return complex


def write_complex(complex: list[float, int, list[int]], filename: str):
    """Write a simplicial complex to a file in the filtration format."""
    with open(filename, 'w') as file:
        for val, dim, verts in complex:
            verts_str = " ".join(str(v) for v in verts)
            file.write(f"{val} {dim} {verts_str}\n")


# Write 2d, 3d, 4d, ..., 10d sphere
# def main():
#     for dim in range(2, 11):
#         d_sphere = build_sphere(dim)
#         write_complex(d_sphere, f"./filtrations/{dim}_sphere")

#     torus = build_torus()
#     write_complex(torus, "./filtrations/torus")

#     filtered_d_sphere_complex = read_filtration("./filtrations/5_sphere")
#     sort_simplices(filtered_d_sphere_complex)
#     print(filtered_d_sphere_complex)
#     boundary_matrix = build_boundary_matrix_sparse(filtered_d_sphere_complex)
#     reduced_matrix, low = reduce_matrix_sparse(boundary_matrix)
#     barcode = compute_barcode(reduced_matrix, low, filtered_d_sphere_complex)
#     write_barcode(barcode, "./output/5_sphere_barcode")
#     plot_barcode(barcode)

#     filtered_torus_complex = read_filtration("./filtrations/torus")
#     sort_simplices(filtered_torus_complex)
#     boundary_matrix = build_boundary_matrix_sparse(filtered_torus_complex)
#     reduced_matrix, low = reduce_matrix_sparse(boundary_matrix)
#     barcode = compute_barcode(reduced_matrix, low, filtered_torus_complex)
#     write_barcode(barcode, "./output/torus_barcode")

# if __name__ == "__main__":
#     main()

# 6. Comparisons TODO:

# For the example of d-dimensional sphere, we get betti numbers 1,0,0,...,0,1 coherent with exercise 2
# For the others, TODO:
# Do we need to check for the whole filtration ?

# 7. Timing the filtrations
import time

def main():
    for file_suffix in ["A", "B", "C", "D"]:
        start_time = time.time()
        filename = f"./filtrations/filtration_{file_suffix}.txt"
        print(f"Processing filtration from file: {filename}")
        
        reading_start_time = time.time()
        filtration = read_filtration(filename)
        sort_simplices(filtration)
        reading_end_time = time.time()
        print(f"Filtration read and sorted in {reading_end_time - reading_start_time:.2f} seconds")

        processing_start_time = time.time()
        boundary_matrix = build_boundary_matrix_sparse(filtration)
        processing_end_time = time.time()
        print(f"Boundary matrix built in {processing_end_time - processing_start_time:.2f} seconds")

        reduction_start_time = time.time()
        reduced_matrix, low = reduce_matrix_sparse(boundary_matrix)
        reduction_end_time = time.time()
        print(f"Boundary matrix reduced in {reduction_end_time - reduction_start_time:.2f} seconds")

        barcode_start_time = time.time()
        barcode = compute_barcode(reduced_matrix, low, filtration)
        barcode_end_time = time.time()
        print(f"Barcode computed in {barcode_end_time - barcode_start_time:.2f} seconds")

        output_filename = f"./output/large_{file_suffix}_barcode"
        write_barcode(barcode, output_filename)

        end_time = time.time()

        print(f"Barcode written to '{output_filename}' in {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    main()

# 8. Interprétations 
# H0: On note que pour les filtrations A, B, C, D il n'y a qu'une barre infinie à la fin des filtrations
# Donc une seule composante connexe
# En fait c'est logique, la offset filtration on commence avec une boule à chaque point, et on connecte petit à petit les boules
# en les faisant grandir
# Donc à la fin on a une seule grosse boule qui englobe tout le nuage de points

# On prend la C
# Je commence à t = 0 juste un peu après (au début c'est juste tous les points qui sont ajoutés petit à petit)
# On a rapidement qu'une seule composante connexe (pareil pour toutes les autres d'ailleurs)
# ça pourrait pas être des sphères arrangées en cercle ? le cercle de diamètre 9 genre et les sphères individuelles de diamètre 1, 
# pendant que ça grandit, ça fait un vide (dimension 2) qui est terminé quand les sphères se rencontrent toutes, 
# et on a un gros cycle jusqu'à 9 (mais que faire du cycle à 1 ?)
# Peut-être qu'il y a un autre plus petit anneau de sphères à l'intérieur qui crée le cycle à 1
# Mais non car on a vite qu'une seule composante connexe
# Non ça marche pas même un premier vide qui nait tout seul puis plein d'autres vides qui naissent et meurent