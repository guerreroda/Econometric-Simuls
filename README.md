# Notes on Regression Discontinuity Design

## Simulations of Discontinuity

This repository simulates an outcome (Y) that is impacted by a Treatment (Treat). When the treatment is sharp, all units above the cutoff (C) are treated. Non-treated are below the cutoff. When the treatment is fuzzy, units are treated according to a probability (Assignment) that is higher above the cutoff and lower below.
The data generation process is influenced by an unobservable correlated with the treatment.

The key identifying assumption of RD is that individuals are comparable around the cutoff, given a bandwidth (B). Therefore, unobservables do not affect the outcome around the cutoff. However, further from the cutoff, unobservables impact the outcome and units are not comparable.

The omitted variable issue above yields biased OLS estimates.
The Sharp RD estimates are unbiased. The FRD design yields an unbiased LATE.

The file RD Simulation.do produces the coefficient in OLS, IV, and RD for a sharp and fuzzy design. The output is a distribution of coefficients as shown below:

<img src="https://github.com/guerreroda/RegressionDiscontinuity/blob/main/RD_Simulation.png" width="50%" height="50%">

The red vertical line is the true coefficient.

## Simulation of Placebo Test

The RD estimator typically uses a cutoff-placebo test as robutsness.

In a sharp design, researchers calculate the correlation between the running variable and the outcome across a series of fake or placebo cutoffs. In the sharp framework, the cutoff separates treated and non-treated observations, such that regressions for any cutoff below the true cutoff include treated observations in the treated group, leading to a smaller or insignificant coefficient than the true cutoff. Placebo cutoffs above the true cutoff include control observations in the treated group, yielding also a smaller coefficient. Therefore, the true cutoff will invariably yield the largest coefficient in the sharp design.

In a fuzzy framework, there are treated units below the cutoff and non-treated units above the cutoff. Because the cutoff changes the probability of treatment, and not the actual assignment to treatment, the coefficients calculated with placebo cutoffs are not systematically smaller than the true cutoff.

This code performs the placebo test in the cutoff for a sharp and fuzzy design. The output for a simulation is below:


<img src="https://github.com/guerreroda/RegressionDiscontinuity/blob/main/RD_Placebo.png" width="50%" height="50%">
