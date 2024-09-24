# Overview
** Unbiased Commercial Premium**
This repository contains the code for our study on optimizing insurance premiums by combining the strengths of Bayesian and commercial premium approaches.

# Proposed Optimal Premium
Unbiasedness, which plays a significant role in insurance science, implies that the expected premium for the current year is a function of explanatory variables.  
Specifically, the expected value of next year's premium, given the explanatory variables, is equal to the assumed distribution function:
$$E[y_{t+1} | X] = \lambda $$ where $ \lambda \equiv exp(X\cdot \beta) $.

In fact,
* Bayesian Premium: This approach efficiently fits previous conditions but can be difficult for policyholders to interpret.
* Commercial Premium: This method is more intuitive, allowing policyholders to easily understand how their claim history affects future premiums  

We propose a method that optimizes the premium by combining the strengths of both approaches. Our solution divides the premium calculation into two distinct components(i.e, $f(\lambda) \cdot g(n_{1:t})$ where $n_{1:t}$ represents the claim history from the 1st to the t-th claim).

<details>
	<summary>Futher study</summary>
  	<div markdown="1">
      Additionally, we show that under certain conditions in a random effects setting, **the Bayesian premium can be equivalent to the commercial premium**, providing a straightforward and effective solution.
  	</div>
</details>

## How to use each file  
The code is divided into two main parts.   
The first part demonstrates the inefficiency of commercial premiums, while the second part compares the performance of various premiums, including Bayesian premiums and two types of commercial premiums.   

The performance results of inefficient commercial premium can be found in Example 6.
And within the Numerical Study section, you can find the performance results of each model in both the Simulation and Real Data Study).

### Performance Measure
* $ \rm{HMSE(y, \hat{y}})}:= \frac{1}{N} \cdot \sum_limit{i=1}^{N}{(\lambda_i \cdot R - \hat{y_i})^2}$
  where $ y \sim F(mean = \lambda \cdot R), R \sim ~ \Pi$  
* $ \rm{MSE(y, \hat{y})}:= \frac{1}{N} \cdot \sum_limit{i=1}^{N}{(y_i - \hat{y_i})^2}$  
* $ DIX(Prem) := \frac{Var(E(Prem|\lambda))}{Var(Prem)}$  


## Usage of File
|name|content|data|
|---|---|---|
|Pois_Gamma_RE.R|define function of all estimations under pois-gamma setting(including optimization, plot etc)||
|NB_Gamma_RE.R|define function of all estimations under NB-gamma setting(including optimization, plot etc)||
|R_example_Biasedness.R|performance of commerical permium for bayesI||
|R_example_Biasedness_bayes2.R|performance of commerical permium for bayesII||
|sim_Pois_Gamma_RE.R|compare performane of each premium(focused on credibility): SPrem, CPrem, GPrem||
|sim_NB_Gamma_RE.R|compare performane of each premium(focused on credibility): SPrem, CPrem, GPrem||
|real_data_Pois_Gamma.R|||
|real_data_NB_Gamma.R|||


