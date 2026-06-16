from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix0 cimport Matrix


cdef class Matrix_generic_dense(Matrix_dense):
    cdef list _entries
    cdef Matrix_generic_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols)
    cdef void _set_to_product_classical(self, Matrix left, Matrix right) except *
