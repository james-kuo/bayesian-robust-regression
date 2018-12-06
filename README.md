I implemented Gibbs Sampler in R to fit Bayesian student-t “Robust” regression with Jeffrey’s prior on a dataset and compared it with Gaussian linear regression with Jeffrey’s prior. Bayesian student-t regression uses the more fat-tailed t-distribution as sampling distribution to adress the outlier problem in Gaussian regression. Artificial noises are added to a dataset under different settings to compare the posterior distributions of parameters of the student-t and Gaussian model. Compared model performances and robustness of posterior distributions.

`Ting_Yuan_Kuo_Final_Project.pdf` is the paper. 

`Analysis.R`, `Analysis2.R`, and `Analysis3.R` are R codes which support the main `Final Project.rmd` markdown file.

`AnalysisMain.RData`, `Corrupt.RData`, `Corrupt2.RData`, and `Corrupt2alt.RData` contains the processed datasets. 

# Main Results
I compared the two models on a dataset with only 1 independent variable. Right panel is when the quadratic term is added in addition to the linear term, and left panel is when only the linear term is present. Lower theta is corresponds to a larger artificial shock added to the dataset. We can see that in the robust model, the regression line is much more stable under any level of shock in the robust regression compared to the gaussian regression.
![Robust vs. Gaussian Model with Exponential noise](https://github.com/james-kuo/bayesian-robust-regression/blob/master/regression_plot.png)

Shown below are the posterior distributions when only linear term is present. We can see that as the level of shock increases, the posterior distribution in the robust model does not change by much.
![Posterior Distributions of Robust vs. Gaussian Model](https://github.com/james-kuo/bayesian-robust-regression/blob/master/posterior_distributions.png)
