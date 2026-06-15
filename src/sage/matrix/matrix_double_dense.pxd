from sage.matrix.matrix_numpy_dense cimport Matrix_numpy_dense
from sage.matrix.matrix0 cimport Matrix


cdef class Matrix_double_dense(Matrix_numpy_dense):
    cdef int _set_to_product_c_impl(self, Matrix left, Matrix right) except -1
