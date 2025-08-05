This repository contains R scripts and associated files relevant to mining, cleaning, and analyzing electronic health record (EHR) data in the AMB_ETL schema as provided by R3.

The primary product of this repo is found in `BacteremiaTreatmentEffects/`. This directory includes scripts, data, and plots that all result in the estimation of average treatment effects (ATE, ATT, ATC, etc.) of alternative antibiotics in the treatment of pathogen-specific bacteremia using EHR data.

Within this directory, `PrimaryScripts/` contains several scripts that:

001 - Build the cohort datasets from raw EHR data (this script requires a connection to the AMB_ETL schema in the Neptune database)

002a - Calculate covariate balancing weights using a variety of methods, including propensity scores calculated using GBMs (though the weights derived from this method do not perform as well as directly computing weights using energy balancing, entropy balancing, CBPS, or stab. balancing methods), 

002b - Assess the balancing effect of the balancing weights (balance diagnostics)

003a - Estimate average treatment effects using the raw data and balancing weights (003b - plot those treatment effects).

* In this section, we estimate the average treatment effect (ATE), average treatment effect on the treated (ATT), average treatment effect on the controls (ATC), average treatment effect in the overlap region (ATO), etc. using both unadjusted models and models adjusted with various weighting strategies to account for confounding.

* Weighting Methods Used:
   * "unw": Unweighted model (no adjustment for covariates - no weights, no other covariates)
   * "uni": Univariate model using only the binary treatment as a predictor
   * "mlt": Multivariate model, including additional covariates that had a post-weight absolute standardized mean difference (SMD) > 0.1

* Balancing/weighting techniques:
   * "energy" - Energy balancing for ATE, ATT, and ATC. See Huling and Mak (2024).
   * "EB" - Entropy balancing for ATT and ATC. See Hainmueller, J. (2012) 'Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies', Political Analysis (Winter 2012) 20 (1): 25–46.
   * "SBW" - Stable balancing weights for ATE, ATT, and ATC. See Chattopadhyay, A., Hase, C. H., and Zubizarreta, J. R. (2020), "Balancing Versus Modeling Approaches to Weighting in Practice," Statistics in Medicine, 39, 3227-3254.
   * "CBPS" - Covariate balancing propensity score estimation for ATE, ATT, and ATC. See Imai, Kosuke and Marc Ratkovic. 2014. “Covariate Balancing Propensity Score.” Journal of the Royal Statistical Society, Series B (Statistical Methodology). http://imai.princeton.edu/research/CBPS.html.
   * Propensity score-based weights (prop scores calculated by GBM in twang package). See Hirano & Imbens (2001): Estimation of causal effects using propensity score weighting., Imbens & Rubin (2015): Causal Inference for Statistics, Social, and Biomedical Sciences. Cambridge University Press.
      * "ATE" - Standard inverse probability of treatment weighting (IPTW). Estimates the average effect of treatment across the full population by reweighting each group to represent the whole population.
      * "ATTw" = Weights control units to resemble the treated group. Treated units get weight 1, controls are weighted by odds of treatment. Focuses on treatment effect among units that actually received treatment.
      * "ATCw" = Weights treated units to resemble the control group. Controls get weight 1, treated are weighted by inverse odds of treatment. Focuses on effect of treatment on the untreated population.
      * "ATO" = Downweights units with extreme propensity scores. Targets the **overlap** population — individuals who could plausibly receive either treatment. Improves stability and reduces extrapolation. See Li, Morgan, & Zaslavsky (2018): Targeted weighting for overlap population. JASA.
      * "ATS" = Balances groups to match marginal treatment assignment probabilities. Used to construct synthetic populations with controlled treatment group proportions. Li & Greene (2013): Efficient estimation for average treatment effects using the propensity score. Journal of Econometrics.
   * *Note that for some methods, no solution is found to the weight optimization problem.*

* Treatment variable note: since we are comparing two antibiotic treatments, we have to arbitrarily assign one antibiotic as "treatment" and the other as "control" since that is how all of these methods view the problem. For each cohort, "treatment" was chosen as the less common antibiotic (e.g., daptomycin for MRSA and *E. faecalis*).

* Outcomes:
   * Binary 30-day mortality
   * Length of stay > median. Note: If a patient died before discharge, their length of stay was imputed as 30 days.
   
* Model Output:
   * Y-axis in plots shows the (weighted) percentage point difference in the outcome between treatment groups (e.g., treated vs. control).
   * Estimates are computed using the avg_comparisons() function from the marginaleffects R package, which provides marginal contrasts for generalized linear models.
   * Each estimate is accompanied by 95% confidence intervals.

* Interpretation:
   * Comparing unw, uni, and mlt variants helps evaluate how adjustment for covariates and balancing weights changes the estimated treatment effect.
   * Multiple weighting schemes are compared side-by-side to assess the robustness and sensitivity of treatment effect estimates to different adjustment methods.
   * Crucially, different weighting methods are more or less **conservative**. In particular, energy balancing leads to near exact balancing as measured by covariate mean (see the plots in 002b_BalanceDiagnostics), but comes at the cost of reduced effective sample size. This tradeoff exists across the board.

004 - [Still a work in progress] - Use causal forests to identify subgroups within each cohort that have significantly different treatment effects. Previous attempts at using causal decision trees yielded interesting, but unstable results. Causal forests provide more stability and might be a useful tool for identifying clinically-meaningful subpopulations in a data-driven way.

Note: no PHI data is stored in this repository. If anything is found that might threaten data governance, please alert me @ `sam.blechman@gmail.com`