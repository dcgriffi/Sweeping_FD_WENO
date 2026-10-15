# Sweeping positivity-preserving FD WENO for the Euler equations

Research code accompanying

> D. Chloe Griffin and Chi-Wang Shu,
> *A sweeping positivity-preserving high-order finite difference WENO scheme for Euler equations*,
> Research in the Mathematical Sciences **13**, 38 (2026).
> https://doi.org/10.1007/s40687-026-00621-2

```bibtex
@article{GriffinShu2026Sweeping,
  author  = {Griffin, D. Chloe and Shu, Chi-Wang},
  title   = {A sweeping positivity-preserving high-order finite difference {WENO} scheme for {E}uler equations},
  journal = {Research in the Mathematical Sciences},
  volume  = {13},
  pages   = {38},
  year    = {2026},
  doi     = {10.1007/s40687-026-00621-2}
}
```

## Scope

This repository contains the research code used to produce the numerical results in
the paper. It is intended as a reference implementation of the method, not as
production software.

- Some parameter values may differ slightly from those stated in the paper due to
  later exploration.
- The code is provided as is and is not actively maintained.

## Motivation

High-order schemes for the compressible Euler equations can produce negative density
or pressure near strong shocks, rarefactions and near-vacuum states, which leads to
non-physical results or blow-up. For finite difference WENO schemes, existing
positivity-preserving techniques either impose strict CFL restrictions (scaling
limiters) or modify the numerical fluxes (flux limiters).

The sweeping procedure is instead a **post-processing step** applied to the solution
after each Runge–Kutta stage, without changing the underlying scheme. It extends the
scalar sweeping technique of Liu, Cheng and Shu [1] from density to the **nonlinear
pressure function**, and more generally to any concave function of the conserved
variables, using the scaling limiter of Zhang and Shu [2]. The procedure:

- is **conservative** and **positivity-preserving** for density and pressure;
- keeps the **high-order accuracy** of the underlying scheme (fifth order here);
- applies to finite difference, finite volume and discontinuous Galerkin schemes;
- is cheap in practice: one or two forward/backward sweeps per Runge–Kutta stage were
  sufficient in all tests, including high-CFL runs at the instability threshold.

## Method in brief

Each solver uses a fifth-order finite difference WENO spatial discretization
(characteristic-wise, Lax–Friedrichs flux splitting) with third-order SSP Runge–Kutta
time stepping. After each RK stage, if any point has negative density or pressure, a
conservative **sweeping** post-processing step is applied: a density sweep, followed by
the pressure sweep with the nodal adjustment described in the paper (forward and
backward sweeps; in 2D, alternating column-by-column and row-by-row "snake" orderings).

Switches near the top of each `main.f90`:

| Parameter | Meaning |
|---|---|
| `use_limiter = 2` | sweeping post-processing (used for all paper results) |
| `use_limiter = 1` | Zhang–Shu positivity-preserving flux limiter, for comparison |
| `swp_num` | maximum number of forward/backward sweeps per RK stage |

At the end of a run, the solvers print the sweeping statistics reported in the paper
(max/average/total sweeps, and max/average/total adjustments with negative and positive
neighbors).

## Examples

| Folder | Paper | Language / parallelism | Setup |
|---|---|---|---|
| `ex1_double_rarefaction_1d/` | Example 1, Fig. 1 | MATLAB | 200 points, T = 0.3 |
| `ex2_sedov_1d/` | Example 2, Fig. 2 | Fortran, serial | 800 points on [−2, 2], T = 0.001 |
| `ex3_vortex_accuracy_2d/` | Example 3, Fig. 3, Table 1 | Fortran, serial | 45²–360² points, T = 0.01 |
| `ex4_sedov_2d/` | Example 4, Fig. 4 | Fortran + MPI (32 ranks) | 640² points on [0, 1.3]², T = 1 |
| `ex5_mach2000_jet/` | Example 5, Fig. 5, Table 2 | Fortran, serial | T = 0.001 |
| `ex6_shock_diffraction/` | Example 6, Fig. 6 | Fortran + MPI (32 ranks) | 1040 × 880 points, T = 2.3 |
| `ex7_detonation/` | Example 7, Fig. 7 | Fortran + MPI (32 ranks) | 400 × 400 points, T = 0.6 |

Tables 3–5 collect the sweeping statistics printed by the runs of all seven examples.

## Building and running

Requirements: `gfortran`; an MPI Fortran compiler (`mpif90`) for examples 4, 6 and 7;
MATLAB for example 1 and the MATLAB plotting scripts. All parameters (grid size, final
time, CFL, ...) are compile-time constants in `main.f90`.

Serial Fortran examples (2, 3, 5):

```bash
make
./a
```

Output files (`den.dat`, `pre.dat`, ...) are written to the working directory.

MPI examples (4, 6, 7) require exactly `cut_x * cut_y = 32` ranks and take the output
directory as their first argument:

```bash
make
mkdir -p out
mpirun -np 32 ./a ./out
```

On a SLURM cluster, `slurm_script.sh` does the same (create `./output` for the job logs
first; results go under `$OUTPUT_ROOT`, default `./runs`). The `module load` lines are
specific to the cluster the runs were made on and will need adjusting.
`ex5_mach2000_jet/slurm_script.sh` runs the serial code as a batch job.

Example 1 (MATLAB): run `FD_WENO_Application.m` (problem 5 is the double rarefaction).
It writes `den.dat`, `pre.dat` and `vel.dat`.

### Reproducing the tables

- **Table 1**: run `ex3_vortex_accuracy_2d` with `nx = 45, 90, 180, 360`. Set
  `new_table = .true.` for the first run (creates `err_table.dat`) and `.false.` for
  the following runs (appends). `plot_err_table.m` prints the table. Each run reads
  the previous resolution's errors from `fort.7` (to compute the order) and then
  overwrites it, so create a placeholder before the first run:
  `echo '  0.10000E+01  0.10000E+01' > fort.7`.
- **Table 2**: the δ comparison was made by changing the cap `0.25` in the nodal
  adjustment of `ex5_mach2000_jet/main.f90` by hand and re-running.
- **Tables 3–5**: printed to standard output at the end of each run.

## Plotting

- `*.lay` files are Tecplot layouts that read the `.dat` files in the working directory.
- `ex2_sedov_1d/plot_data.m`, `ex3_vortex_accuracy_2d/plot_data.m`: MATLAB plots.
  The ex3 script saves figures into `results_2D-vortex-accuracy/`.
- `ex4_sedov_2d/plot_radial.m`: radial projection for Fig. 4. It reads `den_mat.dat`
  (the `den.dat` output without its two Tecplot header lines) and the exact solution
  `sedov.dat`.

### Exact Sedov solutions

The exact Sedov blast-wave solutions used as reference curves in Figs. 2 and 4 are
**not included**. Some plotting files expect them (`density.dat`, `pressure.dat`,
`velocity.dat` for example 2; `sedov.dat` for example 4). They can be generated with a
publicly available Sedov solver, for example the codes accompanying

> J. R. Kamm and F. X. Timmes, *On efficient generation of numerically robust Sedov
> solutions*, Los Alamos National Laboratory report LA-UR-07-2849 (2007),

distributed by F. X. Timmes at https://cococubed.com/research_pages/sedov.shtml.

## References

1. Y. Liu, Y. Cheng, C.-W. Shu, *A simple bound-preserving sweeping technique for
   conservative numerical approximations*, J. Sci. Comput. **73** (2017) 1028–1071.
   https://doi.org/10.1007/s10915-017-0395-x
2. X. Zhang, C.-W. Shu, *On maximum-principle-satisfying high order schemes for scalar
   conservation laws*, J. Comput. Phys. **229** (2010) 3091–3120.
   https://doi.org/10.1016/j.jcp.2009.12.030

## Acknowledgments

The Fortran solvers build on finite difference WENO codes written by Xiangxiong Zhang
(Brown University, 2011), including his parallel (MPI) third-order Runge–Kutta +
fifth-order finite difference WENO solver for the 2D Euler equations, which the MPI
examples extend with the sweeping procedure. We thank him for making these codes
available.

The research of the first author is partially supported by the NSF Graduate Research
Fellowship Program. The research of the second author is partially supported by NSF
grant DMS-2309249 and DOE grant DE-NA0004265.

## License

See [LICENSE](LICENSE).
