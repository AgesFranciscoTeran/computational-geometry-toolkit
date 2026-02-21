# Computational Geometry Toolkit

Implementation of core computational geometry algorithms with emphasis on robustness, correctness, and empirical evaluation.

---

## 1. Geometric Toolkit

Reusable Python module `geometria_toolkit.py` implementing:

- `orient2d(p, q, r)`  
  Robust orientation predicate using adaptive or fallback precision methods.

- `convex_hull(points)`  
  Optimal convex hull algorithm (Monotone Chain / Graham Scan).

- `segment_intersection(segments)`  
  Plane sweep algorithm (Bentley–Ottmann approach).

- Auxiliary geometric utilities.

### Design Goals

- Clean modular implementation  
- Numerical robustness  
- Edge-case handling (collinearity, degeneracy)  
- Unit testing per function (pytest-based)

---

## 2. Smallest Enclosing Circle (SEC)

Implementation of the Minimum Enclosing Circle problem using:

### 2.1 Randomized Welzl Algorithm

Recursive randomized algorithm with boundary set constraint:

- Base case: trivial circle from ≤ 3 boundary points  
- Robust point-in-circle test  
- Degenerate case handling  

### 2.2 Deterministic Alternative

Comparison implementation:

- Convex hull pre-filtering  
- Exhaustive evaluation of pairs and triples  
- Selection of minimal enclosing circle  

### Experimental Comparison

Execution time comparison between:

- Welzl algorithm
- Deterministic hull-based approach

Datasets tested:

- 10
- 100
- 1,000
- 10,000 points

---

## 3. Applied Case Study — Antenna Placement in Quito

Application of the Minimum Enclosing Circle to determine optimal antenna placement covering urban regions in Quito.

Tasks performed:

- Compute minimum enclosing circle
- Visualize:
  - Input points
  - Convex hull
  - Enclosing circle
- Estimate coverage efficiency:
  - Compare circle area vs convex hull area
- Report:
  - Optimal center (geographic coordinates)
  - Radius (km)

![Antenna Coverage Example](report/Figure_1.png)

---

## Technical Stack

- Python 3.10+
- numpy
- matplotlib
- pytest

---

## Repository Structure

```
src/        → Core algorithms
tests/      → Unit tests
report/     → Technical documentation
```

---

## Run Tests

```bash
pytest
```

## Benchmark Summary

Empirical comparison between the randomized Welzl algorithm and the deterministic hull-based approach shows near-linear behavior for Welzl on large datasets (n ≥ 1000), with significant speedup over exhaustive methods.
