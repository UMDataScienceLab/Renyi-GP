# An Alternative Gaussian Process Objective Based on the R\'enyi Divergence

## Introduction

* We introduce an alternative closed-form objective function $\alpha$-ELBO for parameter estimation in the Gaussian process ($\mathcal{GP}$) based on the R\'enyi $\alpha$-divergence. 
* We use a decreasing temperature parameter $\alpha$ to iteratively deform the objective function during optimization. Ultimately, our objective function  converges to the exact likelihood function of $\mathcal{GP}$. At early stages of optimization, $\alpha$-ELBO can be viewed as a regularizer that smoothens out some unwanted critical points. At late stages, $\alpha$-ELBO recovers the exact likelihood function that guides the optimizer to solutions that best explain the observed data. 
* Theoretically, we derive an upper bound of the R\'enyi divergence under the proposed objective and derive convergence rates for a class of smooth and non-smooth kernels. 
* Case studies on a wide range of real-life engineering applications demonstrate that our proposed objective is a practical alternative that offers improved prediction performance over several state of the art inference techniques.

## Prerequisite

* [R](https://www.r-project.org/)

## Code

* main_1.R: the main file that is used to learn our model
* logL.R: return the log-likelihood function
* branin.R, griewank.R, grlee12.r: the simulation functions
