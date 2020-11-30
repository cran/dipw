# Description

Estimation of the average treatment effect when controlling for
    high-dimensional confounders using debiased inverse propensity score
    weighting (DIPW). DIPW relies on the propensity score following a sparse
    logistic regression model, but the regression curves are not required to be
    estimable. Despite this, our package also allows the users to estimate 
    the regression curves and take the estimated curves as input to our
    methods. Details of the methodology can be found in Yuhao Wang and
    Rajen D. Shah (2020) "Debiased Inverse Propensity Score Weighting for
    Estimation of Average Treatment Effects with High-Dimensional Confounders"
    [arXiv link](https://arxiv.org/abs/2011.08661). The package relies on the optimisation
    software [`MOSEK`](https://www.mosek.com/) which must be installed separately;
    see the documentation for `Rmosek`.

# Usage instruction

Once installed, please use `?dipw.ate` and `?dipw.mean` to check the user manual.

# Download

You can download this package via cran, for example using the R command `install.packages("dipw")` in your R console.
