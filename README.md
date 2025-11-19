# T_Diet
T-DIET - Developing novel statistical methods for the analysis of longitudinal dietary patterns and their association with health outcomes
# Multilevel Hidden Markov Models (MHMM) - Comprehensive Mathematical Framework

## Table of Contents
1. [Introduction](#introduction)
2. [Mathematical Foundation](#mathematical-foundation)
3. [Standard Hidden Markov Models](#standard-hidden-markov-models)
4. [Multilevel Hidden Markov Models](#multilevel-hidden-markov-models)
5. [Notation Following Zucchini et al. (2016)](#notation-following-zucchini-et-al-2016)
6. [Model Specification](#model-specification)
7. [Parameter Estimation](#parameter-estimation)
8. [Inference and Decoding](#inference-and-decoding)
9. [Model Selection and Diagnostics](#model-selection-and-diagnostics)
10. [Applications to Nutrition Data](#applications-to-nutrition-data)
11. [Implementation Details](#implementation-details)
12. [Advanced Extensions](#advanced-extensions)
13. [Computational Considerations](#computational-considerations)
14. [References](#references)

---

## Introduction

Hidden Markov Models (HMMs) have become a powerful statistical tool for modeling sequential data with underlying latent states. The extension to **Multilevel Hidden Markov Models (MHMMs)** allows for the analysis of hierarchical or clustered time series data, where observations are nested within higher-level units such as individuals, households, or geographical regions.

This document provides a comprehensive mathematical treatment of MHMMs, following the rigorous notation and framework established in Zucchini, MacDonald, and Langrock's seminal work *"Hidden Markov Models for Time Series: An Introduction Using R"* (2nd edition, 2016).

### Motivation for Multilevel HMMs
---
Traditional HMMs assume independence across observation sequences, which is often violated in practice when:
- Multiple individuals are observed over time
- Repeated measurements are taken from the same experimental units
- Spatial or temporal clustering exists in the data
- Subject-specific heterogeneity needs to be accounted for

MHMMs address these challenges by incorporating random effects and hierarchical structure into the HMM framework.

---

## Mathematical Foundation

### Probability Spaces and Stochastic Processes

Let $(\Omega, \mathcal{F}, \mathbb{P})$ be a probability space. A **stochastic process** $\{X_t\}_{t \in \mathbb{T}}$ is a collection of random variables defined on this space, where $\mathbb{T}$ is typically $\mathbb{N}$ or $\mathbb{Z}^+$ for discrete time.

### Markov Chains

A **Markov chain** $\{S_t\}_{t=1}^T$ of order $m$ with state space $\mathcal{S} = \{1, 2, \ldots, M\}$ satisfies:

$$\mathbb{P}(S_t = s_t | S_{t-1} = s_{t-1}, \ldots, S_1 = s_1) = \mathbb{P}(S_t = s_t | S_{t-1} = s_{t-1}, \ldots, S_{t-m} = s_{t-m})$$

For a **first-order** Markov chain (most common in HMMs):

$$\mathbb{P}(S_t = j | S_{t-1} = i) = \Gamma_{ij}$$

where $\boldsymbol{\Gamma} = [\Gamma_{ij}]_{i,j=1}^M$ is the **transition probability matrix** with $\sum_{j=1}^M \Gamma_{ij} = 1$ for all $i$.

---

## Standard Hidden Markov Models

### Model Definition

A standard HMM consists of:
1. A hidden Markov chain $\{S_t\}_{t=1}^T$ with state space $\mathcal{S}$
2. An observable process $\{Y_t\}_{t=1}^T$ conditional on $\{S_t\}$

The **joint distribution** factorizes as:

$$p(y_{1:T}, s_{1:T}) = \pi_{s_1} \prod_{t=2}^T \Gamma_{s_{t-1},s_t} \prod_{t=1}^T p(y_t | s_t)$$

where:
- $\boldsymbol{\pi} = (\pi_1, \ldots, \pi_M)$ is the initial state distribution
- $\Gamma_{ij}$ are transition probabilities
- $p(y_t | s_t)$ are **emission probabilities** or **conditional distributions**

### Emission Distributions

Common choices for $p(y_t | s_t)$ include:
- **Gaussian**: $Y_t | S_t = j \sim \mathcal{N}(\mu_j, \sigma_j^2)$
- **Poisson**: $Y_t | S_t = j \sim \text{Poisson}(\lambda_j)$
- **Binomial**: $Y_t | S_t = j \sim \text{Binomial}(n, p_j)$
- **Multinomial**: $Y_t | S_t = j \sim \text{Multinomial}(n, \boldsymbol{p}_j)$

---

## Multilevel Hidden Markov Models

### Hierarchical Structure

Consider $K$ level-2 units (e.g., individuals), each with $T_k$ observations at level 1. Let:

- $S^{(k)}_t$ be the hidden state for unit $k$ at time $t$
- $Y^{(k)}_t$ be the observation for unit $k$ at time $t$
- $k = 1, 2, \ldots, K$ index the level-2 units
- $t = 1, 2, \ldots, T_k$ index time within each unit

### Random Effects Approach

In MHMMs, we introduce random effects to capture between-unit heterogeneity:

#### Random Transition Probabilities

Let $\boldsymbol{\Gamma}^{(k)}$ be the transition matrix for unit $k$:

$$\Gamma^{(k)}_{ij} = \frac{\exp(\eta_{ij} + b^{(k)}_{ij})}{\sum_{l=1}^M \exp(\eta_{il} + b^{(k)}_{il})}$$

where:
- $\boldsymbol{\eta} = [\eta_{ij}]$ are fixed effects (population-level transition parameters)
- $\boldsymbol{b}^{(k)} = [b^{(k)}_{ij}]$ are random effects for unit $k$
- $\boldsymbol{b}^{(k)} \sim \mathcal{N}(\boldsymbol{0}, \boldsymbol{\Sigma}_b)$

#### Random Emission Parameters

For Gaussian emissions with unit-specific means:

$$Y^{(k)}_t | S^{(k)}_t = j \sim \mathcal{N}(\mu_j + u^{(k)}_j, \sigma_j^2)$$

where $\boldsymbol{u}^{(k)} = (u^{(k)}_1, \ldots, u^{(k)}_M) \sim \mathcal{N}(\boldsymbol{0}, \boldsymbol{\Sigma}_u)$

### Joint Distribution

The joint distribution for all units is:

$$p(\boldsymbol{y}, \boldsymbol{s}, \boldsymbol{b}, \boldsymbol{u}) = \prod_{k=1}^K \left[ \pi_{s^{(k)}_1} \prod_{t=2}^{T_k} \Gamma^{(k)}_{s^{(k)}_{t-1}, s^{(k)}_t} \prod_{t=1}^{T_k} p(y^{(k)}_t | s^{(k)}_t, \boldsymbol{u}^{(k)}) \right] \times p(\boldsymbol{b}^{(k)}) p(\boldsymbol{u}^{(k)})$$

---

## Notation Following Zucchini et al. (2016)

### Basic Notation

| Symbol | Description |
|--------|-------------|
| $T$ | Length of observation sequence |
| $M$ | Number of hidden states |
| $\mathcal{S} = \{1, \ldots, M\}$ | State space |
| $S_t$ | Hidden state at time $t$ |
| $Y_t$ | Observation at time $t$ |
| $\boldsymbol{\Gamma}$ | Transition probability matrix |
| $\boldsymbol{\pi}$ | Initial state distribution |
| $\boldsymbol{\theta}$ | Emission parameters |
| $\boldsymbol{\delta}$ | Stationary distribution |

### Multilevel Extensions

| Symbol | Description |
|--------|-------------|
| $K$ | Number of level-2 units |
| $T_k$ | Observations for unit $k$ |
| $S^{(k)}_t$ | Hidden state for unit $k$ at time $t$ |
| $Y^{(k)}_t$ | Observation for unit $k$ at time $t$ |
| $\boldsymbol{\Gamma}^{(k)}$ | Transition matrix for unit $k$ |
| $\boldsymbol{b}^{(k)}$ | Random effects for unit $k$ |
| $\boldsymbol{\Sigma}_b$ | Covariance matrix of random effects |

### Matrix Notation

Let $\boldsymbol{P}_t^{(k)}$ be the $M \times M$ diagonal matrix of emission probabilities for unit $k$ at time $t$:

$$\boldsymbol{P}_t^{(k)} = \text{diag}(p(y^{(k)}_t | S^{(k)}_t = 1), \ldots, p(y^{(k)}_t | S^{(k)}_t = M))$$

The **forward probabilities** are computed as:

$$\boldsymbol{\alpha}_t^{(k)} = \boldsymbol{P}_t^{(k)} \boldsymbol{\Gamma}^{(k)T} \boldsymbol{\alpha}_{t-1}^{(k)}$$

with initialization $\boldsymbol{\alpha}_1^{(k)} = \boldsymbol{P}_1^{(k)} \boldsymbol{\pi}$.

---

## Model Specification

### Two-Level MHMM

#### Level 1 (Observation Model)

For unit $k$ at time $t$:

$$Y^{(k)}_t | S^{(k)}_t = j, \boldsymbol{u}^{(k)} \sim f(y | \boldsymbol{\theta}_j + \boldsymbol{u}^{(k)}_j)$$

#### Level 2 (State Process Model)

$$S^{(k)}_t | S^{(k)}_{t-1} = i, \boldsymbol{b}^{(k)} \sim \text{Categorical}(\boldsymbol{\Gamma}^{(k)}_{i\cdot})$$

where $\boldsymbol{\Gamma}^{(k)}_{i\cdot} = g(\boldsymbol{\Gamma}_i, \boldsymbol{b}^{(k)}_i)$

#### Random Effects Distribution

$$\begin{bmatrix} \boldsymbol{b}^{(k)} \\ \boldsymbol{u}^{(k)} \end{bmatrix} \sim \mathcal{N}\left(\boldsymbol{0}, \begin{bmatrix} \boldsymbol{\Sigma}_b & \boldsymbol{\Sigma}_{bu} \\ \boldsymbol{\Sigma}_{ub} & \boldsymbol{\Sigma}_u \end{bmatrix}\right)$$

### Three-Level MHMM

For studies with additional nesting (e.g., measurements within days within individuals):

$$Y^{(l)}_{kt} | S^{(l)}_{kt} = j, \boldsymbol{u}^{(l)}_{kt} \sim f(y | \boldsymbol{\theta}_j + \boldsymbol{u}^{(l)}_{kt})$$

where $l = 1, \ldots, L$ level-3 units, $k = 1, \ldots, K_l$ level-2 units within level-3 unit $l$.

---

## Parameter Estimation

### Maximum Likelihood Estimation

The **likelihood function** for MHMMs involves integrating out random effects:

$$\mathcal{L}(\boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma}) = \prod_{k=1}^K \int \left[ \sum_{s_{1:T_k}} \pi_{s_1} \prod_{t=2}^{T_k} \Gamma^{(k)}_{s_{t-1},s_t} \prod_{t=1}^{T_k} p(y^{(k)}_t | s_t, \boldsymbol{u}^{(k)}) \right] p(\boldsymbol{b}^{(k)}, \boldsymbol{u}^{(k)}) d\boldsymbol{b}^{(k)} d\boldsymbol{u}^{(k)}$$

This integral is typically intractable, requiring approximation methods.

### EM Algorithm

The **Expectation-Maximization** algorithm treats the hidden states and random effects as missing data:

#### E-step: Compute Expected Sufficient Statistics

$$Q(\boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma} | \boldsymbol{\theta}^{(old)}, \boldsymbol{\Gamma}^{(old)}, \boldsymbol{\Sigma}^{(old)}) = \mathbb{E}_{\boldsymbol{s},\boldsymbol{b},\boldsymbol{u}|\boldsymbol{y}, \boldsymbol{\theta}^{(old)}}[\log p(\boldsymbol{y}, \boldsymbol{s}, \boldsymbol{b}, \boldsymbol{u} | \boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma})]$$

#### M-step: Maximize Expected Complete-Data Log-Likelihood

Update parameters by solving:

$$\frac{\partial Q}{\partial \boldsymbol{\theta}} = 0, \quad \frac{\partial Q}{\partial \boldsymbol{\Gamma}} = 0, \quad \frac{\partial Q}{\partial \boldsymbol{\Sigma}} = 0$$

### Monte Carlo EM

When the E-step is intractable, use Monte Carlo integration:

1. Sample $\boldsymbol{s}^{(k,m)}$, $\boldsymbol{b}^{(k,m)}$, $\boldsymbol{u}^{(k,m)}$ from $p(\boldsymbol{s}, \boldsymbol{b}, \boldsymbol{u} | \boldsymbol{y}, \boldsymbol{\theta}^{(old)})$
2. Approximate $Q$ with Monte Carlo average
3. Perform M-step with approximated $Q$

### Direct Maximization via Numerical Integration

Use adaptive quadrature or Laplace approximation:

$$\log \mathcal{L} \approx \sum_{k=1}^K \log \left[ \sum_{s_{1:T_k}} \pi_{s_1} \prod_{t=2}^{T_k} \Gamma_{s_{t-1},s_t} \prod_{t=1}^{T_k} p(y^{(k)}_t | s_t, \hat{\boldsymbol{u}}^{(k)}) \right] - \frac{1}{2} \log|\boldsymbol{\Sigma}| + \text{const}$$

where $\hat{\boldsymbol{u}}^{(k)}$ is the mode of $p(\boldsymbol{u}^{(k)} | \boldsymbol{y}^{(k)}, \boldsymbol{s}^{(k)})$.

---

## Inference and Decoding

### State Decoding

#### Viterbi Algorithm for MHMMs

For each unit $k$, find the most likely state sequence:

$$\boldsymbol{s}^{(k)*}_{1:T_k} = \arg\max_{\boldsymbol{s}^{(k)}_{1:T_k}} p(\boldsymbol{s}^{(k)}_{1:T_k} | \boldsymbol{y}^{(k)}_{1:T_k}, \boldsymbol{\theta}, \boldsymbol{\Gamma})$$

The algorithm proceeds with unit-specific transition matrices $\boldsymbol{\Gamma}^{(k)}$.

#### Local Decoding

Compute posterior state probabilities:

$$\gamma^{(k)}_{tj} = \mathbb{P}(S^{(k)}_t = j | \boldsymbol{y}^{(k)}_{1:T_k}) = \frac{\alpha^{(k)}_{tj} \beta^{(k)}_{tj}}{\sum_{l=1}^M \alpha^{(k)}_{tl} \beta^{(k)}_{tl}}$$

where $\alpha^{(k)}_{tj}$ and $\beta^{(k)}_{tj}$ are forward and backward probabilities for unit $k$.

### Random Effects Estimation

#### Empirical Bayes

Posterior distribution of random effects:

$$p(\boldsymbol{b}^{(k)}, \boldsymbol{u}^{(k)} | \boldsymbol{y}^{(k)}, \boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma}) \propto p(\boldsymbol{y}^{(k)} | \boldsymbol{b}^{(k)}, \boldsymbol{u}^{(k)}, \boldsymbol{\theta}, \boldsymbol{\Gamma}) p(\boldsymbol{b}^{(k)}, \boldsymbol{u}^{(k)} | \boldsymbol{\Sigma})$$

Posterior modes (empirical Bayes estimates):

$$(\hat{\boldsymbol{b}}^{(k)}, \hat{\boldsymbol{u}}^{(k)}) = \arg\max_{\boldsymbol{b},\boldsymbol{u}} p(\boldsymbol{b}, \boldsymbol{u} | \boldsymbol{y}^{(k)}, \boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma})$$

---

## Model Selection and Diagnostics

### Information Criteria

#### Akaike Information Criterion (AIC)

$$\text{AIC} = -2\log \mathcal{L} + 2p$$

where $p$ is the number of parameters:
- Fixed effects: $M(M-1) + (M-1) + M \times \dim(\boldsymbol{\theta})$
- Random effects covariance: $\dim(\boldsymbol{\Sigma})$

#### Bayesian Information Criterion (BIC)

$$\text{BIC} = -2\log \mathcal{L} + p \log(n)$$

where $n = \sum_{k=1}^K T_k$ is the total number of observations.

### Cross-Validation

#### Leave-One-Unit-Out Cross-Validation

For each unit $k$:
1. Fit model on all units except $k$
2. Compute predictive log-likelihood for unit $k$
3. Average across all units

### Goodness-of-Fit Tests

#### Residual Analysis

For Gaussian emissions:

$$r^{(k)}_t = \frac{y^{(k)}_t - \hat{\mu}^{(k)}_{s_t}}{\hat{\sigma}_{s_t}}$$

where $\hat{\mu}^{(k)}_{s_t} = \hat{\mu}_{s_t} + \hat{u}^{(k)}_{s_t}$

#### Pseudo-Residuals

Transform to standard normal under the fitted model:

$$\tilde{r}^{(k)}_t = \Phi^{-1}(F_{\hat{\theta}_{s_t}}(y^{(k)}_t))$$

where $F_{\hat{\theta}_{s_t}}$ is the fitted CDF and $\Phi^{-1}$ is the standard normal quantile function.

---

## Applications to Nutrition Data

### Data Structure

Nutrition studies often have hierarchical structure:
- **Level 1**: Repeated dietary assessments (24-hour recalls, food diaries)
- **Level 2**: Individuals (participants)
- **Level 3**: Households or geographic regions (optional)

### Model Specification for Nutrition Data

#### Food Group Percentages as Observations

Let $Y^{(k)}_{t,g}$ be the percentage energy from food group $g$ for individual $k$ at time $t$:

$$\boldsymbol{Y}^{(k)}_t = (Y^{(k)}_{t,1}, \ldots, Y^{(k)}_{t,G})$$

where $G = 18$ food groups.

#### Multivariate Gaussian Emissions

$$\boldsymbol{Y}^{(k)}_t | S^{(k)}_t = j, \boldsymbol{U}^{(k)} \sim \mathcal{N}_G(\boldsymbol{\mu}_j + \boldsymbol{U}^{(k)}_j, \boldsymbol{\Sigma}_j)$$

with constraint $\sum_{g=1}^G Y^{(k)}_{t,g} = 100\%$.

#### Dirichlet Emissions (Alternative)

$$\boldsymbol{Y}^{(k)}_t | S^{(k)}_t = j, \boldsymbol{U}^{(k)} \sim \text{Dirichlet}(\boldsymbol{\alpha}_j + \boldsymbol{U}^{(k)}_j)$$

### Interpretation of Hidden States

Hidden states represent **dietary patterns**:
- **State 1**: "Health-conscious" pattern (high fruits, vegetables, whole grains)
- **State 2**: "Western" pattern (high processed foods, sugary beverages)
- **State 3**: "Mediterranean" pattern (high olive oil, fish, nuts)

### Time-Varying Covariates

Incorporate seasonality or intervention effects:

$$\text{logit}(\Gamma^{(k)}_{ij}) = \eta_{ij} + b^{(k)}_{ij} + \gamma_{ij} \cdot \text{Season}_t + \delta_{ij} \cdot \text{Intervention}_t$$

---

## Implementation Details

### Computational Algorithms

#### Forward-Backward Algorithm with Random Effects

For each unit $k$:

1. **Forward Pass**:
   $$\alpha^{(k)}_t(s) = p(y^{(k)}_1, \ldots, y^{(k)}_t, S^{(k)}_t = s)$$
   $$\alpha^{(k)}_t = \boldsymbol{P}^{(k)}_t \boldsymbol{\Gamma}^{(k)T} \boldsymbol{\alpha}^{(k)}_{t-1}$$

2. **Backward Pass**:
   $$\beta^{(k)}_t(s) = p(y^{(k)}_{t+1}, \ldots, y^{(k)}_{T_k} | S^{(k)}_t = s)$$
   $$\beta^{(k)}_t = \boldsymbol{\Gamma}^{(k)} \boldsymbol{P}^{(k)}_{t+1} \boldsymbol{\beta}^{(k)}_{t+1}$$

3. **Posterior Probabilities**:
   $$\gamma^{(k)}_t(s) = \frac{\alpha^{(k)}_t(s) \beta^{(k)}_t(s)}{\sum_{l=1}^M \alpha^{(k)}_t(l) \beta^{(k)}_t(l)}$$

#### Numerical Integration Methods

1. **Gauss-Hermite Quadrature**:
   $$\int f(x) \phi(x) dx \approx \sum_{i=1}^q w_i f(x_i)$$
   where $x_i$ and $w_i$ are quadrature points and weights.

2. **Laplace Approximation**:
   $$\int \exp(f(x)) dx \approx \exp(f(\hat{x})) \sqrt{\frac{2\pi}{-f''(\hat{x})}}$$
   where $\hat{x}$ maximizes $f(x)$.

### Software Implementation

#### Template Model Builder (TMB)

TMB provides automatic differentiation for efficient optimization:

```cpp
template<class Type>
Type objective_function<Type>::operator()()
{
  // Data
  DATA_VECTOR(y);           // Observations
  DATA_MATRIX(X);           // Design matrix
  DATA_INTEGER(n_states);   // Number of states
  DATA_INTEGER(n_subjects); // Number of subjects
  
  // Parameters
  PARAMETER_VECTOR(logit_gamma);  // Transition probabilities
  PARAMETER_VECTOR(mu);           // Emission means
  PARAMETER_VECTOR(log_sigma);    // Emission SDs
  PARAMETER_VECTOR(b);            // Random effects
  PARAMETER_MATRIX(Sigma_b);      // Random effects covariance
  
  // Transformations
  vector<Type> gamma = invlogit(logit_gamma);
  vector<Type> sigma = exp(log_sigma);
  
  // Log-likelihood computation
  Type nll = 0;
  // ... forward-backward algorithm implementation
  
  return nll;
}
```

#### R Package Integration

```r
library(hmmTMB)
library(TMB)

# Compile TMB template
compile("mhmm.cpp")
dyn.load(dynlib("mhmm"))

# Parameter initialization
params <- list(
  logit_gamma = rep(0, n_states * (n_states - 1)),
  mu = rep(0, n_states * n_food_groups),
  log_sigma = rep(0, n_states * n_food_groups),
  b = matrix(0, n_subjects, n_random_effects),
  Sigma_b = diag(n_random_effects)
)

# Data list
data <- list(
  y = food_group_percentages,
  X = design_matrix,
  n_states = 3,
  n_subjects = n_individuals
)

# Optimization
obj <- MakeADFun(data, params, DLL = "mhmm")
opt <- nlminb(obj$par, obj$fn, obj$gr)
```

---

## Advanced Extensions

### Time-Varying Transition Probabilities

#### Non-Homogeneous MHMM

Allow transition probabilities to vary with time or covariates:

$$\Gamma^{(k)}_{ij}(t) = \frac{\exp(\eta_{ij} + b^{(k)}_{ij} + \boldsymbol{\gamma}_{ij}^T \boldsymbol{z}^{(k)}_t)}{\sum_{l=1}^M \exp(\eta_{il} + b^{(k)}_{il} + \boldsymbol{\gamma}_{il}^T \boldsymbol{z}^{(k)}_t)}$$

where $\boldsymbol{z}^{(k)}_t$ are time-varying covariates.

### Semi-Parametric Extensions

#### Mixture of Experts Model

$$p(y^{(k)}_t | S^{(k)}_t = j) = \sum_{m=1}^M \pi_{jm} f_m(y^{(k)}_t | \boldsymbol{\theta}_{jm})$$

where $f_m$ are different parametric families (Gaussian, t, Gamma).

#### Bayesian Nonparametrics

Dirichlet Process Prior for emission distributions:

$$G_j \sim DP(\alpha_j, G_0)$$
$$p(y | S_t = j) \sim G_j$$

### Multivariate Extensions

#### Factor-Analytic Structure

For high-dimensional observations:

$$\boldsymbol{Y}^{(k)}_t | S^{(k)}_t = j \sim \mathcal{N}_p(\boldsymbol{\mu}_j + \boldsymbol{\Lambda}_j \boldsymbol{f}^{(k)}_t, \boldsymbol{\Psi}_j)$$

where $\boldsymbol{f}^{(k)}_t \sim \mathcal{N}_q(\boldsymbol{0}, \boldsymbol{I})$ are latent factors.

#### Copula-Based Dependence

Model dependence across food groups using copulas:

$$C(F_1(y_1), \ldots, F_G(y_G) | S_t = j)$$

where $F_g$ are marginal distributions and $C$ is a copula function.

---

## Computational Considerations

### Scalability

#### Computational Complexity

- **Standard HMM**: $O(T M^2)$ per forward-backward pass
- **MHMM**: $O(K T M^2)$ plus integration cost
- **Integration**: $O(q^d)$ where $q$ is quadrature points and $d$ is random effects dimension

#### Memory Requirements

- Store forward/backward matrices: $O(K T M)$
- Random effects: $O(K d)$
- Covariance matrices: $O(d^2)$

### Optimization Strategies

#### Variational Inference

Approximate posterior with factorized distribution:

$$q(\boldsymbol{s}, \boldsymbol{b}, \boldsymbol{u}) \approx q(\boldsymbol{s}) q(\boldsymbol{b}) q(\boldsymbol{u})$$

Minimize KL divergence:

$$\text{KL}(q || p) = \int q \log \frac{q}{p}$$

#### Stochastic Optimization

For large datasets, use stochastic gradient descent:

$$\boldsymbol{\theta}^{(t+1)} = \boldsymbol{\theta}^{(t)} - \eta_t \nabla \mathcal{L}_k(\boldsymbol{\theta}^{(t)})$$

where $\mathcal{L}_k$ is likelihood for a minibatch of units.

### Parallel Computing

#### Embarrassingly Parallel

Units are conditionally independent given parameters:

- Parallelize forward-backward across units
- Parallelize Monte Carlo integration
- Use GPU for matrix operations

#### Message Passing Interface (MPI)

```cpp
// MPI implementation for distributed computing
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

// Distribute units across processes
int units_per_process = K / size;
int start_unit = rank * units_per_process;
int end_unit = (rank + 1) * units_per_process;

// Compute local likelihood
Type local_loglik = compute_loglik(start_unit, end_unit);

// Reduce across processes
Type global_loglik;
MPI_Reduce(&local_loglik, &global_loglik, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
```

---

## References

### Core References

1. **Zucchini, W., MacDonald, I. L., & Langrock, R. (2016).** *Hidden Markov Models for Time Series: An Introduction Using R* (2nd ed.). CRC Press.

2. **Bartolucci, F., Farcomeni, A., & Simonacci, V. (2013).** *Longitudinal Data Analysis Using Hidden Markov Models*. Springer.

3. **Cappé, O., Moulines, E., & Rydén, T. (2005).** *Inference in Hidden Markov Models*. Springer.

### Methodological Papers

4. **Altman, R. M. (2007).** Mixed hidden Markov models: an extension of the hidden Markov model to the longitudinal data setting. *Journal of the American Statistical Association*, 102(477), 201-210.

5. **Maruotti, A. (2011).** Mixed hidden Markov models for longitudinal data: a review. *Statistical Methods in Medical Research*, 20(5), 447-465.

6. **Rydén, T. (2008).** EM versus Markov chain Monte Carlo for estimation of hidden Markov models: A computational perspective. *Brazilian Journal of Probability and Statistics*, 22(2), 140-159.

### Applications in Nutrition

7. **Zhang, S., et al. (2020).** Dietary pattern analysis using hidden Markov models: A novel approach to understanding eating behavior trajectories. *Nutrition Journal*, 19(1), 1-12.

8. **Jiang, Y., et al. (2018).** Multilevel hidden Markov models for analyzing dietary intake patterns in longitudinal studies. *American Journal of Clinical Nutrition*, 108(3), 567-576.

### Computational Methods

9. **Kristensen, K., et al. (2016).** TMB: Automatic differentiation and Laplace approximation. *Journal of Statistical Software*, 70(5), 1-21.

10. **Bürkner, P. C. (2017).** brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80(1), 1-28.

---

## Mathematical Appendices

### Appendix A: Derivation of Forward-Backward Equations

#### Forward Recursion

Starting with $\alpha_1(i) = \pi_i p(y_1 | S_1 = i)$, for $t = 2, \ldots, T$:

$$\alpha_t(j) = p(y_1, \ldots, y_t, S_t = j)$$
$$= \sum_{i=1}^M p(y_1, \ldots, y_{t-1}, S_{t-1} = i, S_t = j, y_t)$$
$$= \sum_{i=1}^M p(y_1, \ldots, y_{t-1}, S_{t-1} = i) \cdot p(S_t = j | S_{t-1} = i) \cdot p(y_t | S_t = j)$$
$$= \sum_{i=1}^M \alpha_{t-1}(i) \Gamma_{ij} p(y_t | S_t = j)$$

#### Backward Recursion

Starting with $\beta_T(i) = 1$, for $t = T-1, \ldots, 1$:

$$\beta_t(i) = p(y_{t+1}, \ldots, y_T | S_t = i)$$
$$= \sum_{j=1}^M p(S_{t+1} = j | S_t = i) \cdot p(y_{t+1}, \ldots, y_T | S_{t+1} = j)$$
$$= \sum_{j=1}^M \Gamma_{ij} p(y_{t+1} | S_{t+1} = j) \beta_{t+1}(j)$$

### Appendix B: EM Algorithm Derivation

#### Complete-Data Log-Likelihood

$$\log p(\boldsymbol{y}, \boldsymbol{s}, \boldsymbol{b}, \boldsymbol{u} | \boldsymbol{\theta}, \boldsymbol{\Gamma}, \boldsymbol{\Sigma}) = \sum_{k=1}^K \left[ \log \pi_{s^{(k)}_1} + \sum_{t=2}^{T_k} \log \Gamma^{(k)}_{s^{(k)}_{t-1}, s^{(k)}_t} + \sum_{t=1}^{T_k} \log p(y^{(k)}_t | s^{(k)}_t, \boldsymbol{u}^{(k)}) + \log p(\boldsymbol{b}^{(k)}, \boldsymbol{u}^{(k)} | \boldsymbol{\Sigma}) \right]$$

#### Expected Sufficient Statistics

Define:
- $\xi^{(k)}_{tij} = \mathbb{P}(S^{(k)}_{t-1} = i, S^{(k)}_t = j | \boldsymbol{y}^{(k)})$
- $\gamma^{(k)}_{ti} = \mathbb{P}(S^{(k)}_t = i | \boldsymbol{y}^{(k)})$

Then:
$$Q = \sum_{k=1}^K \left[ \sum_{i=1}^M \gamma^{(k)}_{1i} \log \pi_i + \sum_{t=2}^{T_k} \sum_{i=1}^M \sum_{j=1}^M \xi^{(k)}_{tij} \log \Gamma_{ij} + \sum_{t=1}^{T_k} \sum_{i=1}^M \gamma^{(k)}_{ti} \log p(y^{(k)}_t | i, \boldsymbol{u}^{(k)}) \right] + \text{random effects terms}$$

### Appendix C: Laplace Approximation

For integral $\int \exp(f(x)) dx$:

1. Find mode: $\hat{x} = \arg\max_x f(x)$
2. Taylor expansion around $\hat{x}$:
   $$f(x) \approx f(\hat{x}) - \frac{1}{2}(x - \hat{x})^2 |f''(\hat{x})|$$
3. Approximate integral:
   $$\int \exp(f(x)) dx \approx \exp(f(\hat{x})) \int \exp\left(-\frac{1}{2}(x - \hat{x})^2 |f''(\hat{x})|\right) dx$$
   $$= \exp(f(\hat{x})) \sqrt{\frac{2\pi}{|f''(\hat{x})|}}$$

---

## Conclusion

Multilevel Hidden Markov Models provide a powerful framework for analyzing hierarchical sequential data, combining the temporal dynamics of HMMs with the flexibility of mixed-effects models. The mathematical foundation, following Zucchini et al.'s notation, offers a rigorous approach to modeling complex real-world phenomena such as dietary patterns, where individual trajectories are nested within population-level processes.

The key advantages of MHMMs include:
- **Flexibility**: Accommodate various data structures and emission distributions
- **Interpretability**: Hidden states represent meaningful latent patterns
- **Efficiency**: Borrow strength across units through random effects
- **Scalability**: Can be extended to high-dimensional and multivariate settings

Future developments in computational methods and software implementations will continue to expand the applicability of these models to increasingly complex datasets in nutrition science and beyond.
