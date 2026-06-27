\name{gse164378_3p_citeseq_hvg}
\alias{gse164378_3p_citeseq_hvg}
\docType{data}
\title{GSE164378 3 Prime CITE-seq RNA HVG and ADT Example}
\description{
A compact paired RNA and ADT example derived from the 3 prime CITE-seq portion
of GSE164378. The object keeps 4,000 matched PBMC cells, 2,000 RNA highly
variable genes, and all 228 measured ADT features.
}
\usage{
data(gse164378_3p_citeseq_hvg)
}
\format{
A list with:
\describe{
\item{rna_counts}{Sparse integer matrix with 2,000 selected RNA genes by 4,000 cells.}
\item{adt_counts}{Sparse integer matrix with 228 ADT features by the same 4,000 cells.}
\item{meta}{Cell metadata, including RNA/ADT QC totals, lane, donor, time, batch, and cell type labels.}
\item{source}{Short source dataset description.}
\item{accession}{GEO accession identifier, \code{"GSE164378"}.}
\item{description}{Short description of the compact paired-modality subset.}
\item{selection}{List describing the cell and RNA feature selection.}
\item{original_dim}{Dimensions of the source 3 prime RNA and ADT matrices before subsetting.}
\item{cell_totals}{RNA and ADT column sums for the selected cells.}
\item{feature_totals}{RNA and ADT row sums for the selected features.}
}
}
\source{
Derived from the 3 prime RNA and ADT matrices in GSE164378, the Hao et al.
multimodal PBMC reference dataset. The source 3 prime matrices have 33,538 RNA
features by 161,764 cells and 228 ADT features by the same 161,764 cells.
}
\examples{
data(gse164378_3p_citeseq_hvg)
rna <- gse164378_3p_citeseq_hvg$rna_counts
adt <- gse164378_3p_citeseq_hvg$adt_counts
dim(rna)
dim(adt)
}
\keyword{datasets}
