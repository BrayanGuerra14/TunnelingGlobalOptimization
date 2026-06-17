# TunnelingGlobalOptimization
# Global Optimization via Steepest Descent and Tunneling

This repository contains a MATLAB implementation of a hybrid global optimization algorithm for multimodal functions. The method combines local minimization by steepest descent, tunneling-based escape mechanisms, and random restarts to explore multiple attraction basins and avoid convergence to suboptimal local minima.

Gradients and directional derivatives are computed through custom automatic differentiation classes based on dual numbers (DualN) and extended dual numbers (DualE), eliminating the need for finite differences. The local optimization stage uses a normalized gradient direction with a derivative-based line search, while the tunneling stage searches for points on the level set associated with the current local minimum using the tunneling function.
