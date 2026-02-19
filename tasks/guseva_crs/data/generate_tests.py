import numpy as np
from scipy.sparse import csr_matrix as crs, random 

# values = [1, 2, 3, 6, 5, 4]
# columns = [0, 1, 2, 0, 1, 2]
# row_inds = [0, 3, 6]

# a = crs((values, columns, row_inds), shape=(2, 3))
# b = crs((
#         [9, 8, 6, 7, 5, 4], 
#         [0, 1, 0, 1, 0, 1], 
#         [0, 2, 4, 6]
#     ), shape=(3, 2))
# c = a.dot(b)
# c.sort_indices()
# print(c.data, c.indices, c.indptr)

def generate_matrix(rows, cols, density):
    return random(rows, cols, density=density, format='csr', data_rvs=lambda s: np.random.randint(1, 100, s) / 2)

def generate_tests(params):
    for size_a, size_b, density, name in params:
        a = generate_matrix(*size_a, density[0])
        b = generate_matrix(*size_b, density[1])
        a.sort_indices()
        b.sort_indices()
        c = a.dot(b)
        c.sort_indices()
        drop_test2file(name + '.txt', a, b, c)
        
def write_matrix(file, a):
    file.write(f"{a.nnz} {a.shape[0]} {a.shape[1]}\n")
    [file.write(f"{str(x)} ") for x in a.data]
    file.write('\n')
    [file.write(f"{str(x)} ") for x in a.indices]
    file.write('\n')
    [file.write(f"{str(x)} ") for x in a.indptr]
    file.write('\n')
        
def drop_test2file(filename, a, b, c):
    with open(filename, 'w', encoding='utf-8') as file:
        write_matrix(file, a)
        write_matrix(file, b)
        write_matrix(file, c)
         
         
        
if __name__ == '__main__':
    generate_tests([
        # ((5, 5), (5, 5), (0.2, 1), 'sparse_dense'),
        # ((10, 10), (10, 10), (1, 0.2), 'dense_sparse'),
        # ((15, 15), (15, 15), (0.1, 0.1), 'double_sparse1'),
        # ((13, 13), (13, 13), (0.1, 0.1), 'double_sparse2'),
        # ((23, 23), (23, 23), (0.1, 0.1), 'double_sparse3'),
        # ((31, 31), (31, 31), (0.1, 0.1), 'double_sparse4'),
        ((5e3, 5e3), (5e3, 5e3), (.005, .005), 'perf')
    ]
    )
    
    