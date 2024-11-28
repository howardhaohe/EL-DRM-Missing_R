# EL-DRM-Missing
An R project for empirical likelihood inference on density ratio model with missing data

This is for my PhD thesis.

"algorithms.R" contains several method to treat missing data:

- Calculating propensity score: kernel method; score method
- Imputation: linear imputation (random, deterministic); kernel imputation; score method imputation

"del.R" contains several ways to define a dual empirical likelihood (DEL) functions for DRM 

Idea 4:
\$
\frac{\delta_{kj}}{\hat{\pi}_{kj}} \frac{\partial L_k(\theta,Y_{kj})}{\partial \beta} + \left( 1-\frac{\delta_{kj}}{\hat{\pi}_{kj}} \right) \frac{\partial L_k(\theta,\hat{Y}_{kj})}{\partial \beta}
\$
