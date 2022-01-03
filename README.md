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


## Installation of PROMES_v05.0.0
1. Go to the lastest release at  https://github.com/saduineveld/Promes_toolbox/releases. Download the zip file. Currently the latest version is PROMES_v05.0.0.a.zip at https://github.com/saduineveld/Promes_toolbox/releases/download/v05.0.0/PROMES_v05.0.0.a.zip;

2. Unpack the zip file in your destination folder. This will add the README.md file, and the folders PROMES_v05.0.0 and TOOLS to the destination folder;

3. To use the Promes toolbox the folder PROMES_v05.0.0 and the subfolders grid_subfun and smolyak_subfun need to be on the searchpath in Matlab. You can run examples from the folder PROMES_v05.0.0/Examples.

All support documents and the licence are found in the the folder PROMES_v05.0.0/DOCUMENTATION.
