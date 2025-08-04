import numpy as np

'''
Bulk modulus by HK
'''

def get_voigt_bulk_modulus(elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the bulk modulus
    K = (elastic_matrix[0, 0] + elastic_matrix[1, 1] + elastic_matrix[2, 2] + 
         2 * (elastic_matrix[0, 1] + elastic_matrix[1, 2] + elastic_matrix[2, 0])) / 9    
    return K

def get_reuss_bulk_modulus(inv_elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if inv_elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the bulk modulus
    K = 1 / (inv_elastic_matrix[0, 0] + inv_elastic_matrix[1, 1] + inv_elastic_matrix[2, 2] + 
             2 * (inv_elastic_matrix[0, 1] + inv_elastic_matrix[1, 2] + inv_elastic_matrix[2, 0]))
    return K

def get_vrh_bulk_modulus(elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the bulk modulus
    K = (get_voigt_bulk_modulus(elastic_matrix) + get_reuss_bulk_modulus(np.linalg.inv(elastic_matrix))) / 2
    return K

'''
Shear modulus
'''

def get_voigt_shear_modulus(elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the shear modulus
    G = ((elastic_matrix[0, 0] + elastic_matrix[1, 1] + elastic_matrix[2, 2]) - 
         (elastic_matrix[0, 1] + elastic_matrix[1, 2] + elastic_matrix[2, 0]) +
         3*(elastic_matrix[3, 3] + elastic_matrix[4, 4] + elastic_matrix[5, 5])
         ) / 15
    return G

def get_reuss_shear_modulus(inv_elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if inv_elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the shear modulus
    G = 15 / (4*(inv_elastic_matrix[0, 0] + inv_elastic_matrix[1, 1] + inv_elastic_matrix[2, 2]) - 
              4*(inv_elastic_matrix[0, 1] + inv_elastic_matrix[1, 2] + inv_elastic_matrix[2, 0]) +
              3*(inv_elastic_matrix[3, 3] + inv_elastic_matrix[4, 4] + inv_elastic_matrix[5, 5])
              )
    return G

def get_vrh_shear_modulus(elastic_matrix):
    # check if the elastic matrix is a 6x6 numpy matrix
    if elastic_matrix.shape != (6, 6):
        raise ValueError("The elastic matrix must be a 6x6 numpy matrix")
    # get the shear modulus
    G = (get_voigt_shear_modulus(elastic_matrix) + get_reuss_shear_modulus(np.linalg.inv(elastic_matrix))) / 2
    return G

'''

'''


if __name__ == '__main__':
    print('This is a module for bulk and shear modulus calculations.')
    mat = np.array([[22, 15, 15, 0, 0, 0], [15, 22, 15, 0, 0, 0], [15, 15, 22, 0, 0, 0], [0,0,0,13,0,0], [0,0,0,0,13,0], [0,0,0,0,0,13]])
    print(get_voigt_bulk_modulus(mat))
    print(get_reuss_bulk_modulus(np.linalg.inv(mat)))
    print(get_vrh_bulk_modulus(mat))
    print('-------------------------------------')
    print(get_voigt_shear_modulus(mat))
    print(get_reuss_shear_modulus(np.linalg.inv(mat)))
    print(get_vrh_shear_modulus(mat))
