Turbulence flow and heat flux budgets for LES/DNS in an OpenFOAM solver
=====================================================================
[![DOI](https://zenodo.org/badge/368479876.svg)](https://zenodo.org/badge/latestdoi/368479876)

The solver is written in OpenFOAM 6. xx, enabling on-the-fly calculation of terms required to calculate budgets. It also contains a postprocessing tool to calculate the necessary gradients for the field's complete flow and heat flux budgets. Please look at test_case and the following commands for an example of running this solver.


The solver and the PostProcess can be built by running the following command in Budgets_PackageV2/PostProcess directory.

```bash
Allwmake
```



To Run the solver, give the following command in Budgets_PackageV2/Solver directory. 

```bash
Allmake
```


In the test case folder, run the following command to get all the results for different entities. 

```bash
Allmake_timefile
```

The solver is named "budgetBuoyantBoussinesqPimpleFoam" 

Controldict and fvschemes files are there for reference.

Please cite this paper if the code is useful to you.

Kakka, Priyesh, and Kameswararao Anupindi. "Assessment of subgrid-scale models for large-eddy simulation of a planar turbulent wall-jet with heat transfer." International Journal of Heat and Mass Transfer 153 (2020): 119593. https://doi.org/10.1016/j.ijheatmasstransfer.2020.119593

References for the equations can be found in "References" directory. Use "Flow budget DNS " for all the flow budget terms, and "Flow budget LES" for extra budget terms due to LES (extra terms are not included in the postProcess but are included in the solver calculation (named "SGS terms"))

I have additional python scripts to create line plots for these budgets (Please contact me at kakkapriyesh@hotmail.com). 

Please let me know if you have any questions. 
