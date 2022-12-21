Turbulence flow and heat flux budgets for LES in an OpenFOAM solver
=====================================================================

The solver is written in OpenFOAM 6. xx, enabling on-the-fly calculation of terms required to calculate budgets. It also contains a postprocessing tool to calculate the necessary gradients for the field's complete flow and heat flux budgets. Please look at test_case and the following commands for an example of running this solver.


The solver and the PostProcess can be built by running the following commands in the post process directory.

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

The solver is named "budgetBuoyantBoussinesqPimpleFoam" 

Controldict and fvschemes files are given for reference.

Please cite this work in the following paper.

Kakka, Priyesh, and Kameswararao Anupindi. "Assessment of subgrid-scale models for large-eddy simulation of a planar turbulent wall-jet with heat transfer." International Journal of Heat and Mass Transfer 153 (2020): 119593. https://doi.org/10.1016/j.ijheatmasstransfer.2020.119593

I have additional python scripts to create line plots for these budgets (Please contact me at kakkapriyesh@hotmail.com) if you want them. 
Please let me know
if you have any questions. 
