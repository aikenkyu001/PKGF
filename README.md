# Parallel Key Geometric Flow (PKGF)
## Purified 12D Manifold Simulation and Linguistic Flow Analysis

All experimental resources are published in this repository: https://github.com/aikenkyu001/PKGF

### Abstract
This repository implements and validates the **Parallel Key Geometric Flow (PKGF)**, a theoretical framework for simulating semantic transitions in a 12-dimensional tangent bundle $TM$. By modeling narrative structures as geometric flows under contextual warping, PKGF ensures the preservation of logical consistency—represented by the **Parallel Key $K$**—through adjoint holonomy updates and co-differential propulsion.

This project provides a rigorous cross-language validation between **Python** (reference implementation) and **Fortran 90/95** (high-performance implementation), demonstrating the conservation of the determinant of $K$ ($\det(K)$) across non-trivial metrics.

---

### Theoretical Framework

#### 1. Geometric Stage
The manifold is defined in $N=12$ dimensions, orthogonally decomposed into four sectors:
$$TM = T_{Subject}M \oplus T_{Entity}M \oplus T_{Action}M \oplus T_{Context}M$$
The metric tensor $g$ is dynamically warped by the intensity of the **Context** sector (Contextual Warping).

#### 2. The Parallel Key $K$
$K \in \Gamma(\mathrm{End}(TM))$ defines the logical structure. It is preserved via parallel transport:
$$K(t+dt) = H K(t) H^{-1}, \quad H = \exp(\Omega dt)$$
where $\Omega$ is the connection matrix derived from the Levi-Civita connection.

#### 3. Fundamental Equations
- **Co-differential Propulsion**: $\frac{\partial}{\partial t}(KX)^\flat = -\delta F = -\star d \star F$
- **Divergence-free Constraint**: $\operatorname{div}_g (KX) = 0$

---

### File Descriptions

| File | Description |
| :--- | :--- |
| `PKGF_Define_jp.md` | Detailed theoretical definitions and mathematical proofs (Japanese). |
| `pkgf_unified_test_en.py` | Python implementation for English linguistic flow analysis. |
| `pkgf_unified_test_jp.py` | Python implementation for Japanese linguistic flow analysis. |
| `pkgf_unified_test_en.f90` | Fortran 90/95 implementation for English flow verification. |
| `pkgf_unified_test_jp.f90` | Fortran 90/95 implementation for Japanese flow verification. |

---

### Prerequisites

- **Python**: 3.x (no external dependencies)
- **Fortran**: `gfortran` (GCC 15.2.0+ recommended)

---

### Building and Running

#### Python Execution
```bash
# Run Japanese analysis
python3 pkgf_unified_test_jp.py

# Run English analysis
python3 pkgf_unified_test_en.py
```

#### Fortran Compilation and Execution
```bash
# Compile and run Japanese version
gfortran -o pkgf_jp pkgf_unified_test_jp.f90
./pkgf_jp

# Compile and run English version
gfortran -o pkgf_en pkgf_unified_test_en.f90
./pkgf_en
```

---

### Scientific Validation
The simulation results across both Python and Fortran confirm:
1. **Conservation of $\det(K)$**: The logical invariant $\det(K) \approx 1.67668$ is preserved with high precision ($< 10^{-15}$ error) during parallel transport.
2. **Behavioral Consistency**: Narrative spikes (e.g., "awakening" events) are correctly detected as sudden inversions in potential and velocity peaks.
3. **Cross-Language Fidelity**: Numerical results for velocity magnitude $||v||$ and divergence $\operatorname{div}_g(KX)$ match between Python and Fortran, with minor variances attributable to IEEE 754 rounding differences in transcendental functions.

---

### License
This research is conducted as part of advanced AGI/ASI studies. All rights reserved.
