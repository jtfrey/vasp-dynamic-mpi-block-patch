# VASP:  Dynamic MPI_BLOCK patch

VASP decomposes some MPI operations into chunks of a fixed size — the `MPI_BLOCK` preprocessor macro — to avoid latency on some architectures.  The source code uses four static arrays (types integer, real, complex, and logical) of dimension `MPI_BLOCK` to buffer reductions, for example.

This patch modifies the source code to dynamically-allocate the four arrays according to the integer-valued `VASP_MPI_BLOCK` environment variable, with the compiled-in `MPI_BLOCK` acting as the default.

The MPI operations that utilize these four arrays are not reentrant:  only one of the arrays is ever in-use at a time.  One existing preprocessor macro enables the reuse of the complex array as storage for the other three arrays via a Fortran equivalence (all other types are smaller than the complex).  Equivalences are not compatible with Fortran allocatables, but the C compatibility routines `C_F_POINTER()` and `C_LOC()` can be used to associate a pointer to the complex array with the other three arrays (as pointers):

```
      COMPLEX(q), DIMENSION(:), ALLOCATABLE :: ZTMP_m
      REAL(q), DIMENSION(:), POINTER :: DTMP_m
      INTEGER, DIMENSION(:), POINTER :: ITMP_m
      LOGICAL, DIMENSION(:), POINTER :: LTMP_m
         :
      ! Pointer-based equivalence between work arrays:
      CALL C_F_POINTER(C_LOC(ZTMP_m), DTMP_m, (/NDTMP/))
      CALL C_F_POINTER(C_LOC(ZTMP_m), ITMP_m, (/NITMP/))
      CALL C_F_POINTER(C_LOC(ZTMP_m), LTMP_m, (/NLTMP/))
```

The arrays at a size of 8000 elements do not put much pressure on memory usage — 125 KiB for the complex array, 62.5 KiB for the real array, and 31.25 KiB for the integer and logical arrays, for a total of 240 KiB.  The equivalence-style reuse of a single memory range woud drop that to 125 KiB for 8000 elements.  This behavior can be forced at compile time by adding the `-DMPI_WORK_ARRAY_EQUIVALENCE` preprocessor option to the VASP build options.

## Applying the patch

Assuming the appropriate patch file exists at `«path-to-patch-file»`, the patch is applied with the `src` subdirectory as the working directory:

```BASH
$ patch -p1 < «path-to-patch-file»
```

The patch can be reversed from the same working directory by adding the `-R` option to the command:

```BASH
$ patch -p1 -R < «path-to-patch-file»
```

## Setting MPI_BLOCK at runtime

With this patch, the runtime dimension of the work arrays is communicated to the program via the `VASP_MPI_BLOCK` environment variable.  If unset, the compiled-in `MPI_BLOCK` is used:

```
   :
 =============================================================
 Using default MPI_BLOCK size 8000 word(s)
 MPI complex work array allocated using MPI_BLOCK
         (with real/integer/logical equivalences)
 running  256 mpi-ranks, with    1 threads/rank, on    4 nodes
 distrk:  each k-point on   64 cores,    4 groups
 distr:  one band on    4 cores,   16 groups
 vasp.6.4.2 20Jul23 (build Aug 25 2023 13:59:04) complex
   :
```

In this case, the `-DMPI_WORK_ARRAY_EQUIVALENCE` option was used at compile time.  Opting for a runtime dimension of 131072 words per array is accomplished by setting `VASP_MPI_BLOCK` in the runtime environment:

```BASH
VASP_MPI_BLOCK=131072 mpirun vasp_std
```

The runtime output reflects this choice:

```
   :
 =============================================================
 Using MPI_BLOCK size 131072 byte(s) from environment
 MPI complex work array allocated using MPI_BLOCK
         (with real/integer/logical equivalences)
 running  256 mpi-ranks, with    1 threads/rank, on    4 nodes
 distrk:  each k-point on   64 cores,    4 groups
 distr:  one band on    4 cores,   16 groups
 vasp.6.4.2 20Jul23 (build Aug 25 2023 13:59:04) complex
    :
 ```
