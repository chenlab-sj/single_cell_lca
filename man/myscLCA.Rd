% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myscLCA.R
\name{myscLCA}
\alias{myscLCA}
\title{LCA scRNA analysis}
\usage{
myscLCA(datmatrix, normalzeddat = NULL, cor.thresh = 0.5,
  clust.max = 10, trainingSetSize = 1000, datBatch = NULL,
  outlier.filter = F, zerocorrection = 0.25)
}
\arguments{
\item{datmatrix}{data matrix, gene by cell}

\item{clust.max}{maximum number of clusters specified by user, default: 10}

\item{trainingSetSize}{the number of cells used in training set; default: 1000; maximum: the total number of cells in data matrix}

\item{datBatch}{a vector of batch group ID for each cell, default: NULL}

\item{zerocorrection}{correction for zeroes before log transformation, default: 0.25}
}
\value{
an integer vector showing clustering membership
}
\description{
This function loads a gene by cell transcript count matarix (row by column).
}
\examples{
# Load example dataset
data(myscExampleData)

# Both data matrix and true labels of cells are provided
names(myscExampleData)

# With 14,074 genes and 250 cells, 83\% of entries are zero
dim(myscExampleData$datamatrix)

# Three types of cells
table(myscExampleData$truelabel)

# Start my scRNAseq LCA analysis
myclust.res <- myscLCA(myscExampleData$datamatrix)

# Top results are provide in a list
length(myclust.res)

# The result ranked as the best, compared with the true labels:
table(myclust.res[[1]],myscExampleData$truelabel)


}
