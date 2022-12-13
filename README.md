# modelling-training-load
Repository for code used in the article "Assessing the cumulative effect of long-term training load on the risk of injury in team sports". 
Available at https://doi.org/10.1136/bmjsem-2022-001342. Modelling the time-varying, cumulative effect of training load

## Overview of scripts
* `functions-relationships.R` has all the functions for the simulated relationship between training load and injury risk.
* `figures-relationships.R` makes figures that illustrate these simulated relationships.
* `simulation-run.R` Uses the above relations, runs a simulation, and all the methods to be compared.
* `functions-performance.R` includes the functions for the performance measures used to compare the methods in simulation-run.R.
* `simulation-results.R` gathers all the results from simulation-run.R to datasets, and runs numerical performance analysis.
* `figures-performance.R` uses the gathered simulation data to illustrate performance with figures.
* `observed-data-analysis.R` is the only script that is not in any way tied to the rest. 
It is a demonstration of the DLNM approach on handball data. These data are not available, but the script shows how it was done for the results in the paper.
