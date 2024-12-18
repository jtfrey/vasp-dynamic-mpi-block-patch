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

If the Fortran compiler can target Fortran 2003 standards, the `HAVE_FORTRAN_2003` macro should be defined in the `Makefile.inc` used to build VASP.

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
 Using MPI_BLOCK size 131072 word(s) from environment
 MPI complex work array allocated using MPI_BLOCK
         (with real/integer/logical equivalences)
 running  256 mpi-ranks, with    1 threads/rank, on    4 nodes
 distrk:  each k-point on   64 cores,    4 groups
 distr:  one band on    4 cores,   16 groups
 vasp.6.4.2 20Jul23 (build Aug 25 2023 13:59:04) complex
    :
 ```

## Altering the M_sum_master algorithm at runtime

The `wave_mpi` module calls the legacy subroutine, `M_sum_master_<type>()`, where `<type>` is `d` or `z` for real versus complex type vectors.  The subroutine performs a reduction of a large vector across all ranks and broadcasts the result (held by the root rank) to all other ranks.  The internal documentation mistakenly labels this an allgather operation:  an allgather merges disparate ranges of data held by ranks into a larger array on all ranks — an `MPI_Gather()` followed by an `MPI_Bcast()`, in essence.  The action in `wave_mpi` is an `MPI_Reduce()` followed by an `MPI_Bcast()`, which equates with an allreduce:  the vector elements are summed across ranks and the results distributed to all ranks.

The motivation for `M_sum_master_<type>()` appears to be the slower, higher-latency networks and ancient MPI libraries that were in use in the early days of HPC.  Likewise, there was no `MPI_Allreduce()` prior to the MPI 2.1 standard so when the code was originally written the `MPI_Reduce()` followed by `MPI_Bcast()` was the norm:  the introduction of `MPI_Allreduce()` was motivated by the prevalence of that code pattern in software.

The code in this patch introduces algorithmic variats in `wave_mpi` and the underlying `M_sum_master_<type>()` controlled by environment variables.  The MPI root rank checks the variables and distributes the state to all other ranks (just as it does with the selected `MPI_BLOCK` size).

| Variable                               | Default     | Discussion                                                                                                                  |
| -------------------------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------- |
| `VASP_MPI_DISABLE_CHUNKED_SUM_MASTER`  | Off/False/0 | When set to On/True/1, the vector will **not** be broken into chunks for the reduction.                                     |
| `VASP_MPI_ENABLE_ALLREDUCE_SUM_MASTER` | Off/False/0 | When set to On/True/1, an `MPI_Allreduce()` will be used in lieu of `MPI_Reduce()` followed by `MPI_Bcast()`. |
