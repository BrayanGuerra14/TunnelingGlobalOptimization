# TunnelingGlobalOptimization
# Global Optimization via Steepest Descent and Tunneling

This repository contains a MATLAB implementation of a global optimization algorithm that combines **Steepest Descent** and **Tunneling Methods** to locate global minima of nonlinear multimodal functions. The implementation employs **automatic differentiation through dual numbers** to compute directional derivatives required during the optimization process.

The algorithm alternates between two stages:

1. **Local Optimization:** A Steepest Descent method is used to locate a local minimum of the objective function.
2. **Tunneling Search:** A tunneling function is constructed around the current minimum to escape its basin of attraction and search for new promising regions of the objective landscape.

This approach allows the algorithm to move beyond locally optimal solutions and improve the exploration of complex optimization landscapes.

## Features

- Steepest Descent local optimization.
- Tunneling-based global search strategy.
- Automatic differentiation using dual numbers.
- Visualization of optimization trajectories.
- Convergence history analysis.
- Benchmark optimization problems.

## Benchmark Functions

The repository includes several well-known test functions used in global optimization:

- Ackley function
- Beale function
- Himmelblau function
- Levy function
- Sphere function

These examples can be selected directly from the main script.

## Repository Structure

```text
.
├── main.m
├── SteepestDescent2.m
├── SteepestDescent3.m
├── grad.m
├── grad_T.m
├── tunnel_c.m
├── DualE.m
└── benchmark_functions/
```

## Requirements

- MATLAB R2020a or newer (recommended).
- No external toolboxes are required beyond the provided automatic differentiation utilities.

## Usage

1. Select the desired benchmark function in `main.m`.
2. Specify the initial point and optimization parameters.
3. Run the main script.

```matlab
fn = @ackley_n;
X0 = [4, 1.5];

main
```

The program generates:

- Optimization trajectories over the objective surface.
- Initial and final points.
- Convergence plots showing the objective value throughout the optimization process.

## Mathematical Background

The implementation combines deterministic local optimization with tunneling techniques inspired by global optimization methods. The local stage follows the negative gradient direction, while the tunneling stage introduces a singular transformation around the current minimum to facilitate transitions toward different attraction basins.

Directional derivatives are computed using dual numbers, avoiding finite-difference approximations and improving numerical accuracy.

## References

- Castellanos, L. *Tunneling Methods for Global Optimization*.
- Griewank, A., and Walther, A. *Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation*.
- Standard benchmark functions from the global optimization literature.

## License

This project is intended for educational and research purposes.
