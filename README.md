# Single cell LCA

A Robust and Scalable Algorithm to Reveal Subtle Diversity in Large-Scale Single-Cell RNA-Seq Data

## Description
Single-cell RNA sequencing (scRNA-seq) emerges as a powerful tool to characterize cell-to-cell variation and dynamics in a seemingly homogenous population.  Efficient and affordable, scRNA-seq is gaining in popularity in both basic and translational biological research areas. However, significant challenges arise in the analysis of scRNA-seq data, including low signal-to-noise ratio with high data sparsity, rising scalability hurdles with hundreds of thousands of cells, and more. Due to inherent complexities in scRNA-seq data, the performance of currently available algorithms may not always be optimal even for fundamental tasks such as identifying heterogeneous subpopulations in the data. In this study, we developed Latent Cellular Analysis (LCA), a machine learning based analytical pipeline that combines similarity measurement by latent cellular states with a graph-based clustering algorithm. LCA features a dual-space model search for both the optimal number of subpopulations and the informative cellular states distinguishing them. LCA provides heuristic solutions for population number inference, dimension reduction, feature selection and confounding factor removal without explicit gene filtering. LCA has proved to be robust, accurate and powerful by comparison to multiple state-of-the-art computational methods on large-scale real and simulated scRNA-seq data. Importantly, LCA’s ability to learn from representative subsets of the data provides scalability, thereby addressing a significant challenge for growing sample size in scRNA-seq data analysis.

## Publications using LCA for single-cell RNAseq analysis

### Metabolic heterogeneity underlies reciprocal fates of TH17 cell stemness and plasticity

[![Chi2019Nature](https://bitbucket.org/scLCA/single_cell_lca/downloads/chi2019nature.png)](https://www.nature.com/articles/s41586-018-0806-7)

Peer W. F. Karmaus, Xiang Chen, Seon Ah Lim, Andrés A. Herrada, Thanh-Long M. Nguyen, Beisi Xu, Yogesh Dhungana, Sherri Rankin, Wenan Chen, Celeste Rosencrance, Kai Yang, Yiping Fan, Yong Cheng, John Easton, Geoffrey Neale, Peter Vogel & Hongbo Chi 

Nature, 2019

Abstract:

	A defining feature of adaptive immunity is the development of long-lived memory T cells to 
	curtail infection. Recent studies have identified a unique stem-like T-cell subset amongst 
	exhausted CD8-positive T cells in chronic infection but it remains unclear whether 
	CD4-positive T-cell subsets with similar features exist in chronic inflammatory conditions. 
	Amongst helper T cells, TH 17 cells have prominent roles in autoimmunity and tissue 
	inflammation and are characterized by inherent plasticity although how such …

### Metabolic signaling directs the reciprocal lineage decisions of αβ and γδ T cells

[![Chi2018ScienceImmuno](https://bitbucket.org/scLCA/single_cell_lca/downloads/chi2019scienceimmun.png)](https://immunology.sciencemag.org/content/3/25/eaas9818.long)

Kai Yang, Daniel Bastardo Blanco, Xiang Chen, Pradyot Dash, Geoffrey Neale, Celeste Rosencrance, John Easton Wenan Chen, Changde Cheng, Yogesh Dhungana, Anil KC, Walid Awad, Xi-Zhi J. Guo, Paul G. Thomas, and Hongbo Chi

Science Immunology, 2018

Abstract:

    The interaction between extrinsic factors and intrinsic signal strength governs thymocyte 
    development, but the mechanisms linking them remain elusive. We report that mechanistic 
    target of rapamycin complex 1 (mTORC1) couples microenvironmental cues with metabolic 
    programs to orchestrate the reciprocal development of two fundamentally distinct T cell 
    lineages, the αβ and γδ T cells. Developing thymocytes dynamically engage metabolic 
    programs including glycolysis and oxidative phosphorylation, as well as mTORC1 signaling …


## Getting Started

### Prerequisites

We need **R** and several **R** packages:

```
	- R
	- R: devtools
	- R: igraph
	- R: mclust
	- R: RMTstat
```


### Installing

Start **R**, then:
```R
	library(devtools)
	install_bitbucket("scLCA/single_cell_lca")
```


### Running the example data in scLCA

In **R**:


```R
# Load the package:

	library(scLCA)

#Load the example dataset provided in the package:

	data(myscExampleData)

#It includes both the transcritp counts matrix and the true labels of the cells:

	names(myscExampleData)
#---output in R---
#[1] "datamatrix" "truelabel"

#With 14,074 genes and 250 cells, 83% of entries are zero

	dim(myscExampleData$datamatrix)
#---output in R---
#[1] 14074	250

# Three types of cells in the example dataset

	table(myscExampleData$truelabel)
#---output in R---
#	1	2	3
#	94	37	119
 
# Start scLCA analysis

	myclust.res <- myscLCA(myscExampleData$datamatrix)
 
# The top result of clustering, compared with the true labels:

	table(myclust.res[[1]],myscExampleData$truelabel)
#---output in R---
#		1	2	3
#	1	94	0	0
#	2	0	1	119
#	3	0	36	0


```


## Corresponding author

###[Xiang Chen](https://www.stjude.org/directory/c/xiang-chen.html)
St. Jude Children's Research Hospital

----





