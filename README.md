# Promes_toolbox
The Promes toolbox solves DSGE models with projection methods, and uses Matlab.

With v05.0.0 you have the options to approximate the policy function with:
- a spline;
- a complete Chebyshev polynomial;
- Smolyak's algorithm (using code by Rafa Valero);
- a complete polynomial based on monomials.

A simple standard RBC model can be solved in less than 0.05 seconds with each of the basis functions. Computation times do increase strongly with the complexity of the model, and a model with four continuous state variables and two policy variables can solved in a couple of seconds.

To use the toolbox the modeler has to carry out four tasks:
- supply a model file which computes the Euler equation residuals;
- set the interval where the policy function should be approximated;
- supply an initial guess for the policy function;
- select an algorithm.

Based on these inputs the toolbox constructs the appropriate grid, and solves the policy function using the selected algorithm. The toolbox includes a function that evaluates the policy function, taking the state variables as input.

Further details are explained in the manual. The toolbox and manual include several examples on how to use the toolbox.
