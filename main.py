import sys

# This code assumes that the filtration is given in an ASCII file with the following format, where each line represents a simplex sigma and is of the form:
# f(sigma) dim(sigma) v_0 ... v_{dim(sigma)}
# where:
# f(sigma) is the function value of sigma (its "time of appearance" in the filtration),
# dim(sigma) is the dimension of sigma,
# v_0 ... v_{dim(sigma)} are the IDs (integers) of the vertices of sigma.
# For instance, 0.125 2 2 6 4 denotes a 2-simplex (triangle) that appears in the filtration at time 0.125 and whose vertices have IDs 2, 6 and 4. Warning: the function values provided in the file must be compatible with the underlying simplicial complex (function values of simplices must be at least as large as the ones of their faces). Nevertheless, the vertex IDs are arbitrary integers and may not start at 0 nor be continuous.

class Simplex:
    def __init__(self, tokens, index):
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


def read_filtration(filename):
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

def sort_simplices(filtration):
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

# 1. Build boundary matrix
def build_boundary_matrix_sparse(filtration):
    """Build the boundary matrix in sparse format from the filtration.
    Returns a list of columns, where each column is represented as a set of row indices with non-zero entries."""
    boundary_matrix = []
    
    for simplex in filtration:
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
                for i, s in enumerate(filtration):
                    if s.dim == simplex.dim - 1 and s.vert == face:
                        boundary.add(i)
                        break
            boundary_matrix.append(boundary)
    
    return boundary_matrix

# Utils: convert sparse matrix to dense matrix
def sparse_to_dense(matrix, n_rows):
    """Convert sparse matrix (list of sets) to dense matrix (list of lists i.e. list of rows)."""
    dense_matrix = []

    for i in range(n_rows):
        dense_matrix.append([0] * len(matrix))
        for j in range(len(matrix)):
            if i in matrix[j]:
                dense_matrix[i][j] = 1
    
    return dense_matrix

# Utils: log dense matrix
def log_dense_matrix(matrix):
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
def reduce_matrix_sparse(matrix):
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