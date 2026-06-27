\name{gse123025_myeloid_hvg}
\alias{gse123025_myeloid_hvg}
\docType{data}
\title{GSE123025 Myeloid HVG Count Matrix Example}
\description{
A larger example subset derived from GSE123025 for RMT spike-model tutorials.
The object keeps all cells from the source matrix and restricts genes to 2,000
highly variable genes selected by \code{select_hvgs()}.
}
\usage{
data(gse123025_myeloid_hvg)
}
\format{
A list with:
\describe{
\item{counts}{Integer matrix with 2,000 selected genes/features by 1,922 cells.}
\item{source}{Source file name in the canonical lab data directory.}
\item{accession}{GEO accession identifier, \code{"GSE123025"}.}
\item{description}{Short description of the HVG subset.}
\item{selection}{List describing the feature and cell selection.}
\item{original_dim}{Dimensions of the source matrix before HVG selection.}
\item{cell_totals}{Column sums of the HVG count matrix.}
\item{gene_totals}{Row sums of the HVG count matrix.}
}
}
\source{
Derived from \code{data/GSE123025/GSE123025_Single_myeloid_1922_cells_processed_data.csv.gz}
in the canonical lab workspace. The source matrix has 23,429 genes by 1,922
cells and is about 8.0 MB compressed as CSV gzip, or 95 MB uncompressed.
}
\examples{
data(gse123025_myeloid_hvg)
counts <- gse123025_myeloid_hvg$counts
dim(counts)
fit <- scscale_fit(counts, r = 10, mp_max_iter = 3, fit_umi = FALSE)
fit
}
\keyword{datasets}
