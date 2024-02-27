## Rectangular Collocation program using LU Decomposition ##

This program provides a small working example for generating the ground state quantum vibrational energy levels of small molecules using the novel rectangular collocation method presented in:

J. Simmons and T. Carrington, "Computing vibrational spectra using a new collocation method with a pruned basis and more points than functions: avoiding quadrature." The Journal of Chemical Physics 158, 144115 (2023) 



### Notes about the program ###

-requires LAPACK + BLAS through the Intel Math Kernel (-mkl) library and ARPACK 

-simply run "make" to build provided required libraries are installed

-"-openmp" flag can be added to LFLAGS for parallel computation

-for more general calculations be cautious that LU decomposition/inverses are numerically stable as the factorization is completed without pivoting

-dneupd2.f is a modified version of the ARPACK dneupd routine that was (at one point) needed by ARPACK for correct eigenvector extraction

-hp_sort2.f is a modified heapsort routine to sort indices in lexicographic order

-uses params.dat file to define potential surface and output file along with pruning parameters and dimension

-the indgen subroutine is used generate the basis, grid, and intermediate indices. Additional indgenND subroutines may be needed depending on the dimension. The nested do-loop form currently used appears faster than a more general index list creator

-the Kinetic energy operator currently uses normal coordinate form, the associated B'' matrices are in the buildmatsmod file

-RAM may be a concern above 9D, the mapping arrays usort,dsort,ssort are large and can likely be improved. Currently due to development time there will be unreferenced (equal to 0) entries owing to structural convenience

