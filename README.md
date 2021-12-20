# Promes_toolbox
The Promes toolbox solves DSGE models with projection methods, and uses Matlab.

With v05.0.0 you have the options to approximate the policy function with:
- spline (see Matlab's griddedInterpolant);
- Complete Chebyshev polynomial;
- Smolyak's algorithm (using code by Rafa Valero);
- Complete polynomials with monomial basis.


To use the toolbox the modeler has to carry out four tasks:
- supply a model file which computes the Euler equation residuals;
- set the interval where the policy function should be approximated;
- supply an initial guess for the policy function;
- select an algorithm.
