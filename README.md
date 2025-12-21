# Evaluation of influenza A/H3N2 epidemiology in England during the 2025-26 season
This repository contains all of the code and data required to reproduce the analyses presented in the accompanying [manuscript](https://zenodo.org/records/17704679).

The analyses should be run in order using the numbered files in the `scripts` folder. Most R packages used are standard, listed at the top of each script, and can be installed from CRAN. However, the [EpiStrainDynamics](https://github.com/acefa-hubs/EpiStrainDynamics/tree/main) package requires some custom input to solve an issue currently raised [here](https://github.com/acefa-hubs/EpiStrainDynamics/issues/34).

Please download the [EpiStrainDynamics](https://github.com/acefa-hubs/EpiStrainDynamics/tree/main) package locally and replace L236 in `R/post_fit_calc.R`:

`a <- transform_posterior_multi(post, B_true, components$num_path, components$num_days)`

with:

`a <- transform_posterior_single(post, B_true, components$num_days)`

All graphical outputs are saved to the `figures` folder, and tables of results are saved either to the `figures` folder or the `results` folder.

The results from the scenario analyses presented in the main text are saved in the `scenarios` folder and are extracted using the scripts `4.combine_scenario_plots.R` and `5.combine_scenario_tables.R`.


## Seasonal influenza H3N2 scenario analysis.

This repo provides analysis and code used to explore and model seasonal H3N2 influenza outbreaks occurring in England, whilst accounting for UK half-term and Christmas holidays alongside measures of immune susceptibility and transmission. Scenario analysis was performed using an age-stratified SIR model that draws contact frequencies from the POLYMOD study using the [`socialmixr`](<https://github.com/epiforecasts/socialmixr>) R package.

The easiest way to run the model is to visit the Shiny app hosted at:

<https://hay-idd.shinyapps.io/ModelFluUk-H3N2/>

You can also run the app locally by running the app.R file in the `ModelFluUk-H3N2` directory.

The app allows users to modify various parameters and explore the resulting impact on flu timing, epidemic size and peak incidence across age groups.


