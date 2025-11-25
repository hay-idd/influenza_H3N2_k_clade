# Seasonal influenza H3N2 scenario analysis.

This repo provides analysis and code used to explore and model seasonal H3N2 influenza outbreaks occurring in England, whilst accounting for UK half-term and Christmas holidays alongside measures of immune susceptibility and transmission. Scenario analysis was performed using an age-stratified SIR model that draws contact frequencies from the POLYMOD study using the socialmixr R package (<https://github.com/epiforecasts/socialmixr>).

For full details please see our pre-print:

<https://www.medrxiv.org/>

### Performing scenario analysis using the Shiny app.

The easiest way to run the model is to visit the Shiny app hosted at:

<https://hay-idd.shinyapps.io/ModelFluUk-H3N2/>

You can also run the app locally by running the app.R file in the `ModelFluUk-H3N2` directory.

The app allows users to modify various parameters and explore the resulting impact on flu timing, epidemic size and peak incidence across age groups. The parameters which can be modified are explained below:

- **`Basic reproductive number R0`** *(float, default=2)*


