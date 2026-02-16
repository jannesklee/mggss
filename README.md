# MGGSS - Multigrid Gauss-Seidel Solver

MGGSS is a C implementation of a multigrid solver for the 2D Poisson equation using the Gauss-Seidel method. It employs Dirichlet boundary conditions on an equidistant 2D grid.

```
Version : 0.0
Author  : Jannes Klee
Contact : jannes.klee@gmail.com
```

## Mathematical Problem

MGGSS solves the **2D Poisson equation** with Dirichlet boundary conditions:

**∇²u = f(x,y)** on the unit square Ω = [0,1] × [0,1]

With boundary conditions:
- **u(x,y) = 0** on ∂Ω (Dirichlet boundary conditions)

The Poisson equation is discretized using **finite differences** on an equidistant grid with step size h = 1/(n+1), where n is the number of interior grid points in each dimension.

### Discrete Formulation

The Laplacian ∇²u is approximated using the **5-point stencil**:

```
∇²u₁ⱼ ≈ (uᵢ₊₁ⱼ + uᵢ₋₁ⱼ + uᵢⱼ₊₁ + uᵢⱼ₋₁ - 4uᵢⱼ) / h² = fᵢⱼ
```

This leads to the linear system:

```
4uᵢⱼ - uᵢ₊₁ⱼ - uᵢ₋₁ⱼ - uᵢⱼ₊₁ - uᵢⱼ₋₁ = h²fᵢⱼ
```

## Numerical Methods

### Gauss-Seidel Smoother

The solver uses the **Gauss-Seidel method** as a smoother, which iteratively updates the solution:

```
uᵢⱼ⁽ᵏ⁺¹⁾ = 0.25 * (h²fᵢⱼ + uᵢ₊₁ⱼ⁽ᵏ⁾ + uᵢ₋₁ⱼ⁽ᵏ⁺¹⁾ + uᵢⱼ₊₁⁽ᵏ⁾ + uᵢⱼ₋₁⁽ᵏ⁺¹⁾)
```

### Multigrid Method

MGGSS implements a **recursive multigrid V-cycle** with the following components:

1. **Pre-smoothing**: Apply Gauss-Seidel iterations on the fine grid
2. **Restriction**: Compute residuals and restrict to coarser grid using full weighting
3. **Coarse grid correction**: Recursively solve the error equation on coarser grids
4. **Prolongation**: Interpolate corrections back to finer grids using bilinear interpolation
5. **Post-smoothing**: Apply Gauss-Seidel iterations on the fine grid

The multigrid method significantly accelerates convergence compared to single-grid methods by efficiently reducing both high-frequency and low-frequency error components.

### Restriction Operator

The restriction operator uses **full weighting** to transfer residuals from fine to coarse grids:

```
v_c[i/2+j/2*(n_c+2)] = 0.25*(v[i+j*(n+2)] + 0.5*(v[(i+1)+j*(n+2)] + v[(i-1)+j*(n+2)] + v[i+(j+1)*(n+2)] + v[i+(j-1)*(n+2)]) + 0.5*(v[(i+1)+(j+1)*(n+2)] + v[(i-1)+(j-1)*(n+2)]))
```

### Prolongation Operator

The prolongation operator uses **bilinear interpolation** to transfer corrections from coarse to fine grids:

- **Even-even points**: Direct copy from coarse grid
- **Even-odd/odd-even points**: Linear interpolation
- **Odd-odd points**: Bilinear interpolation

## Implementation Details

### Data Structures

- **Grid**: Contains solution array `u`, right-hand side array `v`, grid size `n`, and step size `h`
- **Arrays**: Use row-major ordering with ghost points for boundary handling

### Algorithm Parameters

- **n**: Number of grid points (must be odd for proper coarsening)
- **gamma**: Number of multigrid levels (controls recursion depth)
- **eps0**: Error tolerance for convergence
- **k**: Problem type selector (1 = constant density circle, 2 = Gaussian distribution)

### Error Measurement

The solver uses the **maximum norm** of the residual as the convergence criterion:

```
||r||_∞ = max |rᵢⱼ| over all grid points
```

## Usage

### Compilation

```sh
$ make
```

### Cleanup

```sh
$ make clean
```

### Configuration

Setup and initialization can be modified in `src/init.c`:

- **Grid size**: Change `n` variable (must be odd)
- **Error tolerance**: Modify `eps0` 
- **Problem type**: Set `k` to 1 (circle) or 2 (Gaussian)
- **Multigrid levels**: Adjust `gamma` parameter

### Initial Conditions

The solver supports two built-in right-hand side functions:

1. **Constant density circle** (k=1):
   ```
   f(x,y) = 100 if √((x-0.5)² + (y-0.5)²) < 0.1
   f(x,y) = 0   otherwise
   ```

2. **Gaussian distribution** (k=2):
   ```
   f(x,y) = exp(-1000 * ((x-0.5)² + (y-0.5)²))
   ```

## Output Format

The solver outputs the solution in CSV format:

```
x, y, u(x,y), f(x,y)
```

Where:
- `x, y`: Coordinates in [0,1] × [0,1]
- `u(x,y)`: Computed solution
- `f(x,y)`: Right-hand side function

## Mathematical Background

### Multigrid Efficiency

The multigrid method achieves **O(n)** complexity for solving elliptic PDEs, compared to O(n²) or O(n³) for traditional iterative methods. This is possible because:

1. **Smoothing property**: Gauss-Seidel efficiently reduces high-frequency errors
2. **Approximation property**: Coarse grids provide good approximations to low-frequency errors
3. **Two-grid convergence**: The combination of smoothing and coarse grid correction leads to rapid error reduction

### Convergence Analysis

For the Poisson equation with Dirichlet boundary conditions, the multigrid V-cycle typically achieves:

- **Convergence rate**: ~0.1-0.3 per cycle (depending on γ)
- **Asymptotic complexity**: O(n) operations for n unknowns
- **Memory requirements**: O(n) storage

## License

This software is licensed under the GPLv3. See LICENSE for details.

## References

The implementation follows standard multigrid algorithms as described in:

- Briggs, W.L., Henson, V.E., McCormick, S.F. (2000). "A Multigrid Tutorial"
- Hackbusch, W. (1985). "Multi-Grid Methods and Applications"
- Trottenberg, U., Oosterlee, C.W., Schüller, A. (2001). "Multigrid"
