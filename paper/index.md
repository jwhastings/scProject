# **Methods**

In this section we review our datasets, batch correction methods, batch correction evaluations, 

## **Datasets**

A combination of simulated datasets and real scRNA-seq data were used to evaluate the performance of the batch correction methods.

### **Simulated Data**

> Zappia, Luke, Belinda Phipson, and Alicia Oshlack. "Splatter: simulation of single-cell RNA sequencing data." Genome biology 18.1 (2017): 174.
>
> [Splatter: simulation of single-cell RNA sequencing data][paper_splatter]

Datasets 1, 2, ... were simulated using the `splatter` R package.

#### **Dataset 1**



#### **Dataset 2**



### **Real Data**



## **Batch Correction Methods**

We selected <three??> batch correction methods and compared their performance against <five??> different scRNA-seq datasets. A description of each method is listed below.

### **JIVE**

> Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.
> 
> [Joint and individual variation explained (JIVE) for integrated analysis of multiple data types][paper_JIVE]

The joint and individual variation explained** (JIVE) method seeks to decompose two or more biological datasets into three low-rank approximation components: a joint structure among the datasets, individual structures unique to each distinct dataset, and residual noise. If we let $X_1$, $X_2$, ..., $X_k$ be matrices of dimension $p_i \times n$ containing the original datasets, then we have

$$
\begin{align*}
    X_1 &= J_1 + A_1 + \epsilon_1 \\
    X_2 &= J_2 + A_2 + \epsilon_2 \\
    \vdots \\
    X_k &= J_k + A_k + \epsilon_k 
\end{align*}
$$

where $J_i$ denotes the $i$th joint structure submatrix, $A_i$ denotes the $i$th individual structure matrix, and $\epsilon_i$ are error matrices with independent entries.

JIVE was originally created to decompose any set of related biological data: an example from the paper uses gene expression and miRNA data from a set of 234 Glioblastoma Multiforme tumor cells. Our main interest is how well the method performs when we apply it in the context of scRNA-seq data batch correction. We expect the joint structure to serve as the batch-corrected dataset which can then be used in further downstream analyses.

The JIVE decomposition estimates the joint and individual structures by minimizing the sum of squared error. Given an initial estimate for the joint structure, it finds the individual structures to minimize the sum of squared error. Then, given the new individual structures, it finds a new estimate for the joint structure which minimizes the sum of squared error. This process is repeated until a given threshold for convergence is reached. The low-rank approximations are estimated by one of two different methods: a permutation test rank selection and a BIC rank selection. 

#### **Algorithm Improvements**

> O’Connell, Michael J., and Eric F. Lock. "R. JIVE for exploration of multi-source molecular data." Bioinformatics 32.18 (2016): 2877-2879.
>
> [R.JIVE for exploration of multi-source molecular data][paper_R.JIVE]

The JIVE algorithm was implemented into the `r.jive` R package by Dr. Michael O'Connell from the original MATLAB code. 

### **Seurat**

> Stuart, Tim, et al. "Comprehensive integration of single-cell data." Cell 177.7 (2019): 1888-1902.
>
> [Comprehensive Integration of Single-Cell Data][paper_Seurat]

The Seurat integration method employs a sophisticated strategy which "anchors" diverse datasets together. First, log-normalization is performed on all datasets and expression values are standardized for each gene. A subset of features are selected which exhibit high variance across all datasets. Then an initial dimension reduction method utilizing canonical correlation analysis (CCA) is performed to ensure similarities across datasets are preserved. Canonical correlation vectors (CCV) are then approximated and used to identify K-nearest neighbors (KNN) for each cell within their paired dataset. Mutual nearest neighbors (MNN) are then identified to act as anchors between datasets. These anchors are then filtered, scores, weighted, and used to perform the batch correction. 

## **Batch Correction Evaluation**

We employed three metrics to evaluate the performance of each batch correction methods: visual inspection of t-distributed stochastic neighbor embedding (t-SNE) and uniform manifold approximation and projection (UMAP) dimension reduction plots, k-nearest neighbor batch effect tests (kBET), and average silhouette width (ASW).

For the visual inspections, we expect to see cells from different batches overlapping each other in our plots. This is indicative of well-mixed (i.e., integrated) batches.

### **t-Distributed Stochastic Neighbor Embedding**

> Van der Maaten, Laurens, and Geoffrey Hinton. "Visualizing data using t-SNE." Journal of machine learning research 9.11 (2008).
>
> [Visualizing Data using t-SNE ][paper_t-SNE]

t-SNE is a non-linear dimension reduction technique that aids in visualizing high-dimensional data by assigning each data point a location in a two or three-dimensional map. It aims to preserve as much of the local structure of the original data as possible while also revealing global structure such as clusters. High dimensional Euclidean distances between points are used to create conditional probabilities of one point picking the other as its neighbor. A similar conditional probability is calculated for a low dimensional version of the data. The goal of t-SNE is to find a low dimensional (i.e., two or three dimension) representation that matches the two probabilities as best as possible by minimizing a certain cost function. 

### **Uniform Manifold Approximation and Projection**

> McInnes, Leland, John Healy, and James Melville. "Umap: Uniform manifold approximation and projection for dimension reduction." arXiv preprint arXiv:1802.03426 (2018).
>
> [UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction][paper_UMAP]



### **k-Nearest Neighbor Batch Effect Test**

> Büttner, Maren, et al. "A test metric for assessing single-cell RNA-seq batch correction." Nature methods 16.1 (2019): 43-49.
>
> [A test metric for assessing single-cell RNA-seq batch correction][paper_kBET]
> 
> [Github][github_kBET]



### **Average Silhouette Width**

>
>
>



---

[paper_splatter]: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0
[paper_JIVE]: https://doi.org/10.1214/12-aoas597
[paper_R.JIVE]: https://academic.oup.com/bioinformatics/article/32/18/2877/1744043
[paper_Seurat]: https://doi.org/10.1016/j.cell.2019.05.031
[paper_t-SNE]: https://www.jmlr.org/papers/v9/vandermaaten08a.html
[paper_UMAP]: https://doi.org/10.48550/arXiv.1802.03426
[paper_kBET]: https://www.nature.com/articles/s41592-018-0254-1
[github_kBET]: https://github.com/theislab/kBET
