
Turbulence flow and heat flux budgets for LES in an OpenFOAM solver
=====================================================================

The solver is written in OpenFOAM 6.xx and it enables on the fly calculation of terms required to calculate budgets. It also contains a postprocessing tool to
enable the calculation of gradients requried for calculating complete flow and heatflux budgets in the entire field. Please look at test_case and following commands for an example on how to run this solver.

Please cite this work with the following paper

Kakka, Priyesh, and Kameswararao Anupindi. "Assessment of subgrid-scale models for large-eddy simulation of a planar turbulent wall-jet with heat transfer." International Journal of Heat and Mass Transfer 153 (2020): 119593.


The solver and the PostProcess can be built by running following commands in the postprocess directory.

```bash
Allwmake
```



To Run the solver, give the following command in the solver folder

```bash
Allmake
```


In the test case folder, run the following command to get all the results for different entities. 

```bash
Allmake_timefile
```

Solver is named as "budgetBuoyantBoussinesqPimpleFoam" 

Controldict and fvschemes files are given for reference.

I have additional python scripts to create line plots for these budgets (Please contact me at kakkapriyesh@hotmail.com) if you want it. 
Please let me know
if you have any questions. 








