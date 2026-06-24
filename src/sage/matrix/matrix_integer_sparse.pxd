from sage.libs.gmp.types cimport mpz_t
from sage.modules.vector_integer_sparse cimport mpz_vector
from sage.ext.mod_int cimport mod_int
from sage.matrix.matrix_sparse cimport Matrix_sparse

cdef class Matrix_integer_sparse(Matrix_sparse):
    cdef mpz_vector* _matrix

    cdef int mpz_height(self, mpz_t height) except -1
    cdef _mod_int_c(self, mod_int p)
