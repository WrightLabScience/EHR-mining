This repository contains R scripts and associated files relevant to mining, cleaning, and analyzing electronic health record (EHR) data in the AMB_ETL schema as provided by R3.

The primary product of this repo is found in `BacteremiaTreatmentEffects/`. This directory includes scripts, data, and plots that all result in the estimation of average treatment effects (ATE, ATT, ATC, etc.) of alternative antibiotics in the treatment of pathogen-specific bacteremia using EHR data.

Within this directory, `PrimaryScripts/` contains several scripts that:

001 - Build the cohort datasets from raw EHR data (this script requires a connection to the AMB_ETL schema in the Neptune database)

002a - Calculate covariate balancing weights using a variety of methods, including propensity scores calculated using GBMs (though the weights derived from this method do not perform as well as directly computing weights using energy balancing, entropy balancing, CBPS, or stab. balancing methods), 

002b - Assess the balancing effect of the balancing weights (balance diagnostics)

003a - Estimate average treatment effects using the raw data and balancing weights (003b - plot those treatment effects).

004 - [Still a work in progress] - Use causal forests to identify subgroups within each cohort that have significantly different treatment effects. Previous attempts at using causal decision trees yielded interesting, but unstable results. Causal forests provide more stability and might be a useful tool for identifying clinically-meaningful subpopulations in a data-driven way.

Note: no PHI data is stored in this repository. If anything is found that might threaten data governance, please alert me @ `sam.blechman@gmail.com`