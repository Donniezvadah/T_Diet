Here is a comprehensive, mathematically rigorous **README.md** file designed to accompany the R script. It bridges the gap between the raw code and the theoretical framework of the T-DIET project.

***

# T-DIET Project: Longitudinal Dietary Pattern Analysis via Multi-Level Hidden Markov Models

## 1. Overview
This repository contains the implementation of the **T-DIET analytical framework** for identifying latent dietary patterns from longitudinal, compositional food intake data.

Standard statistical methods fail when applied to dietary data due to two fundamental constraints:
1.  **The Compositional Nature:** Dietary data (percentage energy contributions) reside on the simplex $\mathcal{S}^D$, where variables are mutually dependent (they must sum to 100%).
2.  **Longitudinal Dependency:** Observations within subjects over time are not independent; they exhibit temporal autocorrelation and subject-specific heterogeneity.

This codebase solves these problems by combining **Compositional Data Analysis (CoDA)** with **Multi-Level Hidden Markov Models (ML-HMMs)** using the `hmmTMB` engine.

---

## 2. Mathematical Framework

### 2.1 The Sample Space: The Simplex
Let $\mathbf{x} = [x_1, x_2, \dots, x_D]$ be a vector of energy contributions from $D=18$ food groups. The sample space is the simplex:
$$
\mathcal{S}^D = \left\{ \mathbf{x} \in \mathbb{R}^D \mid x_i > 0, \sum_{i=1}^D x_i = \kappa \right\}
$$
where $\kappa=1$ (or 100%). Standard Euclidean geometry (and thus standard covariance matrices) is invalid on $\mathcal{S}^D$.

### 2.2 The Isometry: ILR Transformation
To enable Gaussian modeling, we map the simplex to real Euclidean space $\mathbb{R}^{D-1}$ using the **Isometric Log-Ratio (ILR)** transformation.
Let $\mathbf{z} \in \mathbb{R}^{D-1}$ be the transformed vector:
$$
\mathbf{z} = \text{ilr}(\mathbf{x}) = \mathbf{V}^T \ln(\mathbf{x})
$$
where $\mathbf{V}$ is a $D \times (D-1)$ orthonormal basis matrix satisfying $\mathbf{V}^T \mathbf{V} = \mathbf{I}_{D-1}$ and $\mathbf{V} \mathbf{V}^T = \mathbf{I}_D - \frac{1}{D}\mathbf{1}_D \mathbf{1}_D^T$.

### 2.3 The Hidden Markov Model
We assume the existence of $K$ latent dietary patterns (Hidden States). Let $S_{it} \in \{1, \dots, K\}$ denote the state of subject $i$ at visit $t$.

**The State Process (Markov Chain):**
The transition between dietary patterns is governed by a transition probability matrix $\boldsymbol{\Gamma}$:
$$
\gamma_{jk} = P(S_{it} = k \mid S_{i,t-1} = j)
$$

### 2.4 Multi-Level Emission Probabilities
The core innovation of the T-DIET framework is the inclusion of random effects to account for inter-subject variability. The observation $\mathbf{z}_{it}$ (the PC-reduced ILR coordinates) is modeled as:

$$
\mathbf{z}_{it} \mid (S_{it}=k, \mathbf{u}_i) \sim \mathcal{N}(\boldsymbol{\mu}_k + \mathbf{u}_i, \boldsymbol{\Sigma}_k)
$$

Where:
*   $\boldsymbol{\mu}_k$ is the population-level mean dietary profile for pattern $k$.
*   $\mathbf{u}_i \sim \mathcal{N}(\mathbf{0}, \boldsymbol{\Omega})$ is the **subject-specific random intercept**. This vector shifts the mean, allowing individual $i$ to have a personalized baseline diet while still belonging to global pattern $k$.
*   $\boldsymbol{\Sigma}_k$ is the state-specific covariance.

### 2.5 Estimation
The marginal likelihood is obtained by integrating out the random effects $\mathbf{u}_i$. As this integral is intractable, the code uses the **Laplace Approximation** via the `TMB` (Template Model Builder) backend:

$$
L(\boldsymbol{\theta}) \approx \prod_{i=1}^N \frac{f(\mathbf{z}_i \mid \mathbf{u}_i^*, \boldsymbol{\theta}) p(\mathbf{u}_i^*)}{|H(\mathbf{u}_i^*)|^{1/2}}
$$

---

## 3. Code Structure

The analysis is contained in a single R script (`TDIET_Analysis.R`) organized into four logical modules:

### Module 1: Data Ingestion & Cleaning
*   **Input:** Raw text stream (Nutritics data format).
*   **Dictionary:** Maps `FG1`...`FG18` to human-readable labels (e.g., "Grains", "Meat").
*   **Sanitization:** Removes rows with missing compositional data (`NA`s) and ensures `ID` and `Visit` are treated as factors.

### Module 2: T-DIET Transformation Pipeline
1.  **Closure:** Enforces the constraint $\sum x_i = 1$.
2.  **Zero Handling:** Implements Bayesian-Multiplicative Replacement (`cmultRepl`) to handle zeros (which define the boundary of the simplex), ensuring $\ln(x_i)$ is defined.
3.  **ILR:** Maps $\mathcal{S}^{18} \to \mathbb{R}^{17}$.
4.  **PCA:** Performs dimensionality reduction on the ILR coordinates to retain the top 3 Principal Components ($PC_1, PC_2, PC_3$) for numerical stability in the HMM.

### Module 3: `hmmTMB` Modeling
*   **Distribution:** Multivariate Normal (diagonal covariance structure for PCs).
*   **Formula:** `mean = ~ s(IDNumber, bs = "re")`. This `mgcv`-style syntax instructs the model to estimate a random intercept for every unique `IDNumber`.
*   **Optimization:** Fits the model using the EM algorithm combined with TMB's automatic differentiation.

### Module 4: Publication-Quality Visualization
Generates vector graphics (.pdf) adhering to strict academic standards:
*   **Colorblind-safe palettes** (Okabe-Ito derivative).
*   **Serif fonts** matching LaTeX standard text.
*   **No chartjunk:** Minimalist themes with high data-ink ratios.

---

## 4. Outputs

The script generates the following PDFs in the working directory:

1.  **`1_Dietary_Pattern_Profiles.pdf`**:
    *   *Mathematical meaning:* The geometric mean of the compositions conditioned on the decoded state: $\exp(E[\ln(\mathbf{x}) \mid S=k])$.
    *   *Usage:* Used to label states (e.g., "Prudent Pattern", "Western Pattern").

2.  **`2_PCA_State_Separation.pdf`**:
    *   *Mathematical meaning:* A projection of $\mathbf{z}_{it}$ onto the first two eigenvectors of the covariance matrix, colored by the Viterbi-decoded state sequence.
    *   *Usage:* Visual validation of cluster separability.

3.  **`3_Dietary_Transitions_Alluvial.pdf`**:
    *   *Mathematical meaning:* Visualizes the flux of subjects between states $S_t$ and $S_{t+1}$ over visits $t=1..4$.
    *   *Usage:* Illustrating longitudinal dietary stability vs. volatility.

4.  **`4_Validation_Fibre_Distribution.pdf`**:
    *   *Usage:* External validation. Compares a nutrient *not* explicitly used in the clustering (Fiber) across the derived states to check biological plausibility.

5.  **`5_Missing_Data_Map.pdf`**:
    *   *Usage:* Exploratory Data Analysis (EDA) to assess the sparsity of the input matrix.

---

## 5. Dependencies

The code relies on the following R packages:

```r
install.packages("hmmTMB")       # HMM with Random Effects
install.packages("compositions") # Aitchison Geometry
install.packages("zCompositions")# Zero Replacement
install.packages("ggalluvial")   # Sankey Diagrams
install.packages("RColorBrewer") # Palettes
```

## 6. References

1.  **Aitchison, J. (1986).** *The Statistical Analysis of Compositional Data*. Chapman and Hall.
2.  **Michelot, T. (2022).** *hmmTMB: Hidden Markov Models with Flexible Covariate Effects in R*. arXiv preprint arXiv:2211.14139.
3.  **Zucchini, W., et al. (2016).** *Hidden Markov Models for Time Series: An Introduction Using R*. CRC Press.