# An Alternative Gaussian Process Objective Based on the R\'enyi Divergence

## Introduction

* We introduce an alternative closed-form objective function $\alpha$-ELBO for parameter estimation in the Gaussian process ($\mathcal{GP}$) based on the R\'enyi $\alpha$-divergence. 
* We use a decreasing temperature parameter $\alpha$ to iteratively deform the objective function during optimization. Ultimately, our objective function  converges to the exact likelihood function of $\mathcal{GP}$. At early stages of optimization, $\alpha$-ELBO can be viewed as a regularizer that smoothens out some unwanted critical points. At late stages, $\alpha$-ELBO recovers the exact likelihood function that guides the optimizer to solutions that best explain the observed data. 
* Theoretically, we derive an upper bound of the R\'enyi divergence under the proposed objective and derive convergence rates for a class of smooth and non-smooth kernels. 
* Case studies on a wide range of real-life engineering applications demonstrate that our proposed objective is a practical alternative that offers improved prediction performance over several state of the art inference techniques.

## Prerequisite

* [R](https://www.r-project.org/)

## Files

* CMAPSSData.zip -- the NASA condition monitoring dataset
* contour_plot_1.R, contour_plot_2.r and high_dimension_plot.R -- generate contour plots in our main paper
* different_length.R -- optimizing the $\alpha$-ELBO$ using high-dimensional input data. Each dimension has one length parameter
* experiment.R -- optimizing the $\alpha$-ELBO$ using high-dimensional input data. All dimensions have the same length parameter
* griewank.R and grlee12.r -- simulation functions
* logL.R -- our objective function
* mBCG_noT.R -- the blackbox matrix-matrix multiplication (BBMM) algorithm
* turbine.R -- code to test our algorithm on the NASA condition monitoring dataset

## Implementation

* Our objective function $\alpha$-ELBO is defined in the file "LogL.R". This file includes the implementation of Eq. (7) in Sec. 4.2 in our main paper. 
* Files "different_length.R", "experiment.R", and "turbine.R" contain implementation of our algorithms using different dataset. We evaluate the gradient of $\alpha$-ELBO numerically using the "nloptr" package. In the optimization process, we take a batch of training data and evaluate the numerical gradient of $\alpha$-ELBO. We then update our model parameter by running gradient method implementated in the "nloptr" package.
* The prediction is done through the BBMM algorithm. We implement a conjugate gradient method (in the file "mBCG_noT.R") to solve a quadratic optimization problem in Sec. 6 and obtain the predictive mean and variance.
