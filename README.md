# 3D Pore Boundary Benchmark — MATLAB

**Research code for a comparative study of seven boundary definitions in three-dimensional pore spaces.**

This repository contains **34 original MATLAB `.m` files** recovered from the project's October–November 2026 research materials. The MATLAB algorithms have been preserved without rewriting them. One archived file, originally named `analyser_tore_128_3_points(1).m`, was renamed to `analyser_tore_128_3_points.m` to match its primary MATLAB function.

> **Status:** Experimental research code archive. The files and primary function names have undergone static checks, but **the full pipeline has not been executed or validated in MATLAB** as part of this archive preparation. The source datasets are not included. The different `main_...` scripts represent historical experimental variants and should not be assumed to produce identical results.

## Research overview

A binary 3D image describes the pore space (`1`, `true`, or `255`) and the solid phase (`0`). Seven boundary point sets are compared:

| Method | Definition |
| --- | --- |
| `I6`, `I18`, `I26` | Centers of **pore voxels** adjacent to the solid phase under 6-, 18-, or 26-neighborhood connectivity |
| `E6`, `E18`, `E26` | Centers of **solid voxels** adjacent to the pore phase under 6-, 18-, or 26-neighborhood connectivity |
| `F` | Centers of faces separating pore and solid voxels |

**Domain boundary convention:** The region outside the input image is considered solid. A virtual solid layer is used during boundary extraction; consequently, some external voxel centers may have coordinate `0` or `dim+1`, and face centers may have coordinate `0.5` or `dim+0.5`.

The main computational workflow is:

**Boundary point extraction → Delaunay tetrahedralization → inadmissible tetrahedron removal → circumsphere computation → TD (deviation/spillover) filtering → TC (overlap) selection → reconstruction and geometric evaluation.**

## Source files

The `.m` files are kept at the repository root because several routines invoke dependencies by their MATLAB function names rather than package-qualified paths.

### Core pipeline

| File | Purpose |
| --- | --- |
| `extraire_7_frontieres.m` | Extracts the seven boundary point sets |
| `analyser_delaunay_7_frontieres.m` | Computes Delaunay tetrahedralizations and tetrahedron-quality statistics |
| `nettoyer_tetra_7_frontieres.m` | Filters tetrahedra based on quality and geometric/interface criteria |
| `calculer_boules_7_frontieres.m` | Computes, deduplicates, and classifies circumspheres |
| `analyser_boules_7_frontieres_etape3.m` | Examines the resulting sphere banks |
| `filtrer_boules_TD_7_frontieres.m` | Applies the TD sphere-deviation threshold |
| `selection_TC_volumique_7_frontieres.m` | Performs volume-based TC overlap selection |
| `couverture_union_avant_TC_7_frontieres.m` | Measures union coverage before TC selection |

### Benchmark drivers and historical variants

- `main_benchmark_7_frontieres_64.m` — 64³ benchmark using a `k=4` subvoxel grid.
- `main_benchmark_7_frontieres_64_geoTD.m` — 64³ variant with geometric TD evaluation.
- `main_benchmark_7_frontieres_512.m` — historical 512³ variant; potentially very memory-intensive.
- `main_benchmark_7_frontieres_512_reprise.m` — resume/recovery driver for the 512³ experiment.
- `main_benchmark_7_frontieres_step1_rapide.m` — `k=1` pipeline with voxel-center TD; includes TD0/5/10, TC0/5/10 and an optional geometric TD0 check.
- `main_benchmark_7_frontieres_step1_sequentiel_audit.m` — sequential pipeline with TD audits.
- `main_benchmark_7_frontieres_step1_cubes_intersectes.m` — alternative volumetric evaluation protocol.
- `main_benchmark_7_frontieres_k1_zero.m`, `main_benchmark_7_frontieres_k1_zero_TD0_5_10.m`, and `main_benchmark_7_frontieres_k2_zero_TD0_5_10_TC0_5_10.m` — additional experimental configurations.

### Geometry checks and evaluation

- `verifier_contenance_geometrique_boules_retenues.m`
- `auditer_TD_volume_boules_retenues.m`
- `auditer_TD_continu_7_frontieres.m`
- `auditer_TD_continu_7_frontieres_v2.m`
- `budget_fixe_7_frontieres.m`
- `lancer_budget_fixe_BD12_DRP372.m`
- `calculer_MSD_RSAE_7_frontieres.m`
- `01_surface_BD12.m`, `02_surface_DRP372.m`, `03_surface_Tore.m`
- `validation_physique_boules_DRP372.m`

### Synthetic torus and topology

- `generer_tore_128_benchmark.m`
- `analyser_tore_128_3_points.m`
- `analyser_topologie_tore_3_configurations.m`
- `analyser_surface_tore_TD0_TC0.m`
- `analyser_surface_tore_TD5_TC5.m`

## Requirements

- A recent MATLAB installation with **Image Processing Toolbox**.
- **Statistics and Machine Learning Toolbox** for scripts that call `knnsearch`.
- **Parallel Computing Toolbox** for variants that use `parfor`. Not all scripts provide a sequential fallback.
- Sufficient RAM and disk space for larger volumes, particularly 256³ and 512³. Global tetrahedralizations and subvoxel grids can be computationally demanding.

No specific MATLAB release has been established as fully validated for this archive.

## Quick start

In MATLAB, change to the repository directory and add it to your path:

```matlab
addpath(pwd);

% Small synthetic test; 1 = pore and 0 = solid.
img = false(16,16,16);
img(4:13,4:13,4:13) = true;

% Test boundary extraction without running the heavy benchmark:
[P, stats] = extraire_7_frontieres(img);
disp(stats);
```

To run the complete **64³ benchmark** with a properly prepared pore mask `img`:

```matlab
assert(isequal(size(img), [64 64 64]));
resultats = main_benchmark_7_frontieres_64(img, 'resultats_delaunay_64');
```

For a 3D volume of another size, using the `k=1` protocol:

```matlab
% Use a new or empty output directory for this entry point.
resultats = main_benchmark_7_frontieres_step1_rapide(img, 'new_benchmark_results');
```

To generate the synthetic 128³ torus **without** running a benchmark:

```matlab
[imgTorus, ~] = generer_tore_128_benchmark(false);
```

To generate it **and** run the corresponding benchmark:

```matlab
[imgTorus, torusResults] = generer_tore_128_benchmark(true);
```

## Important input convention: BD12 and DRP-372

The benchmarking functions expect **`0 = solid`, `1/255 = pore`**. However, some original RAW/MAT volumes encode the **pore as `0`**. In that case, invert the binary representation before running the pipeline:

```matlab
img = uint8(originalVolume == 0) * 255;
```

For DRP-372, `02_surface_DRP372.m` reads the `bin` variable from `374_02_024_256.mat` and treats `bin == 0` as pore. For BD12, `01_surface_BD12.m` expects a file named `bd12.mat` with a variable named `img` by default. Adapt filenames, variables, and relative paths to your local environment. The surface-analysis scripts expect previously generated benchmark outputs.

Reference image sizes and pore-voxel counts from the research experiments:

| Dataset | Volume size | Pore voxels |
| --- | ---: | ---: |
| BD12 | 512³ | 23,056,832 |
| DRP-372 | 256³ | 688,309 |
| Synthetic torus | 128³ | 182,384 |

Original image volumes and large generated outputs are not distributed in this repository.

## Reproducibility and limitations

- Detailed algorithms and options are documented in the MATLAB source-file headers. Some implementations use a tetrahedron quality threshold of `q = 10^-3` and a sphere-deduplication tolerance of `10^-7` voxel; these settings are not necessarily universal across all variants.
- **Voxel-based TD0 is not equivalent to exact geometric TD0.** A sphere that contains no solid voxel center may still extend geometrically into the solid phase. This distinction is explicitly addressed in the `step1_rapide` variant.
- Results depend on the benchmark driver and its parameters (`k=1`, `k=2`, `k=4`, continuous audits, etc.). For scientific reproducibility, always report the exact `main_...` script and configuration used.
- Some scripts have relative input paths or expected dataset filenames that must be adjusted to your machine.
- The fixed-budget BD12 loader checks for a variable with `nnz(x) == 23056832`. If your original file uses `0 = pore`, convert the volume or adapt the loader before using it; otherwise, the expected variable may not be detected.
- **No MATLAB execution or independent reproduction of published result tables** was performed when preparing this repository. Validate code paths and outputs before claiming reproducibility.
- Do not commit large `.mat` or `.raw` datasets, reconstructed volumes, or generated figures without first checking distribution rights and provenance.

## Uploading to GitHub

Create a new repository, for example `3d-pore-boundary-benchmark`. Extract the release archive, open a terminal in the `Frontieres3D_GitHub` directory, and run:

```bash
git init
git add .
git commit -m "Initial MATLAB boundary benchmark source"
git branch -M main
git remote add origin https://github.com/<username>/<repository>.git
git push -u origin main
```

Add a `LICENSE` file only after confirming distribution rights for all contributions. Before making the repository public, provide the appropriate author/contributor credits and a definitive bibliographic reference to the associated study, if available.

## Citation

**Citation details to be added after the associated study's bibliographic information and author-approved attribution have been confirmed.**

## License

