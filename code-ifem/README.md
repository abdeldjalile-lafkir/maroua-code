# IFEM Code — 4-Domain Cleanup Report

## Objective

Remove all FEM/domain-5 related code from the original codebase, keeping only the 4 IFEM infinite simplices, while ensuring the code runs correctly and produces valid results.

---

## Problems Encountered & Solutions

### 1. Array Sizing Issue (Root Cause)

- **Problem**: The original uses `nbnodes = 5*...` (6750 nodes for 5 domains). When we downsized to `4*...` (5400 nodes), the code crashed due to buffer overflow.
- **Root cause**: The mesh generation routines (`first_enumeration_whsp`, `categories_nodes_whsp`, `new_number_whsp`) all assume 5 domains internally. Downsizing these arrays corrupts memory.
- **Fix**: Keep `nbnodes = 5*...` = 6750 in `fondements_keltoum.h`, but change `nbtetra = 4*...` = 32000 only (tetrahedra count for 4 IFEM domains). The mesh generation pipeline stays at 5-domain sizes, but tetrahedra are built for only 4 domains.

### 2. Buffer Underflow in `v_rhs(m-1)` (`direct_rhs` line 374)

- **Problem**: Accessing `v_rhs(m-1)` where `m = tetra(k, nb)` can be 1 (the origin node), causing `v_rhs(0)` which is out of bounds.
- **Fix**: Added guard `if (m.ne.1)` before accessing `v_rhs(m-1)`.

### 3. Domain 4 Excluded from Weighting/Error Computation

- **Problem**: Domain 4 (IFEM) was excluded from error computation and RHS generation.
- **Fix**: Changed `dom.lt.4` to `dom.lt.5` in both `weighting` and `example_and_rhs`.

### 4. `hmax(5)` Out-of-Bounds Read

- **Problem**: Original reads `hmax(5)` (domain 5) but only 4 domains exist.
- **Fix**: Changed `ndom=5` to `ndom=4` in `calc_pas` call, removed `hmax(5)` read statement.

### 5. `nbsimp = 5` in `simplexes_keltoum.h`

- **Problem**: Despite removing domain 5, the code needs `nbsimp = 5` to store reference vertices for the 4 IFEM domains + reference big tetrahedron.
- **Fix**: Keep `nbsimp = 5` unchanged.

### 6. Interpolation Crash (`find_big_domain`)

- **Problem**: `find_big_domain` returns garbage when a point is outside all 4 IFEM domains. The faulty `dom.gt.0` guard passed garbage values to `loc_interp_ifem`, crashing the program.
- **Fix**:
  - Initialize `nm = 0` in `find_big_domain`
  - Changed guard to `dom.ge.1 .and. dom.le.4`
  - Skip `find_num_tetr` when no domain found (`s.le.0`)

### 7. "Impossible subdivision of a quadrilateral" Error

- **Root cause**: This was a consequence of the array sizing bug. Corrupted mesh data caused the intersection routine to fail. Fixed automatically when array sizes were corrected.

### 8. Restoring `interpol` and `edge` Subroutines

- Added the missing `interpol` subroutine (from `interpolation.f` in the maroua project), cleaned of domain-5 code.
- Restored the missing `edge` subroutine in `solveur_keltoum.f`.

---

## Final File Structure

| File                         | Role                                                 |
| ---------------------------- | ---------------------------------------------------- |
| `fondements_keltoum.h`       | `nbnodes = 6750`, `nbtetra = 32000`                  |
| `simplexes_keltoum.h`        | `nbsimp = 5` (unchanged)                             |
| `mesh_infinite_keltoum.f`    | Generates 5-domain data but only 4-domain tetrahedra |
| `ordre0_directRhs_keltoum.f` | `m.ne.1` guard for `v_rhs` access                    |
| `interpolation.f`            | `dom.ge.1 .and. dom.le.4` guard                      |
| `infinite_keltoum.f`         | Main program, `ndom=4`, labels changed to `IFEM4`    |
| `calc_erreurs_keltoum_ok.f`  | Error computation for 4 domains only                 |

---

## Results

- **N=20**: 32000 tetrahedra, 6750 nodes (5400 actually used, rest unused)
- **N=6**: 864 tetrahedra, 205 nodes
- Program runs to completion with no errors (Exit: 0)
- All computation phases work: mesh generation → volume check → solution computation → linear solve → error computation
- Error labels: `IFEM` for domains 1-3, `IFEM4` for domain 4

---

## Technical Notes

1. `nbnodes` kept at original 5-domain size (6750) because the numbering/deduplication algorithms assume 5 domains internally. This wastes some memory (1350 unused slots out of 6750) but ensures stability.
2. Tetrahedra are built for only 4 domains (32000), so the linear system only solves for actual IFEM nodes.
3. `find_big_domain` only searches domains 1-4 (IFEM). Points outside these domains are silently skipped during interpolation.

---

## Appendix: Comparison of the Three Codebases

| Property                 | `original/`            | `code/`                       | `code-ifem/`       |
| ------------------------ | ---------------------- | ----------------------------- | ------------------ |
| Number of domains        | 5 (4 IFEM + 1 FEM)     | 5 (4 IFEM + 1 FEM)            | **4 IFEM only**    |
| `nbtetra`                | `5*(3n+N)`             | `4*(3n+N)`                    | `4*(3n+N)`         |
| `nbtetra1`               | —                      | `4*(3n+N)+nt`                 | —                  |
| `nbnodes`                | `5*(...)` = 6750       | `5*(...)` = 6750              | `5*(...)` = 6750   |
| `nbnodes1`               | —                      | 6750 + gm                     | —                  |
| `tetra1` / `xin1`        | —                      | Exists (stores domain 5 data) | —                  |
| Domain 5 source          | `mesh_simp` (analytic) | `construir` + `.dat` files    | **Does not exist** |
| `ndom` in `calc_pas`     | 5                      | 5                             | **4**              |
| `make_tetra` loop        | 5 domains              | 4 + manual domain 5 append    | 4 only             |
| Error label for domain 5 | `MEF`                  | `MEF`                         | —                  |
| Extra files              | —                      | `construir.f`, `.dat` files   | `README.md`        |

### Architecture Evolution

```
original/ : Full 5-domain code (IFEM + FEM), everything generated analytically
    │
    ▼
code/     : Same 5 domains, but domain 5 (FEM) read from external files
            via construir.f, with extra arrays (tetra1, xin1) for extensibility
    │
    ▼
code-ifem/: Domain 5 completely removed — 4 IFEM domains only
            nbnodes=6750 preserved for internal numbering algorithm compatibility
            Tetrahedra built for 4 domains only (nbtetra=32000)
```

---

## Build and Run

```bash
# Build
gfortran -O2 -o code_ifem *.f

# Run (N=20)
./code_ifem

# Run with different N (edit N in fondements_keltoum.h)
sed -i 's/N = 20/N = 6/' fondements_keltoum.h
gfortran -O2 -o code_ifem *.f && ./code_ifem
sed -i 's/N = 6/N = 20/' fondements_keltoum.h

# Debug build
gfortran -g -fcheck=all -fbacktrace -o code_ifem_dbg *.f && ./code_ifem_dbg

# Clean
rm -f code_ifem code_ifem_dbg code_ifem_asan res_form_fin_*.dat essai.dat exact2.dat interpol2.dat fort.* nodes10.dat xi.dat tetra.dat
```
