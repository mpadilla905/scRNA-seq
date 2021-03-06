## ----'quick_intro_01', message = FALSE-------------------------------------------------------------------------------
library('scRNAseq')
library('scater')
library('scran')
library('plotly')


## ----'quick_intro_02', cache = TRUE----------------------------------------------------------------------------------
sce <- scRNAseq::MacoskoRetinaData()

## How big is the data?
pryr::object_size(sce)

## How does it look?
sce


## ----'quick_intro_03', cache = TRUE----------------------------------------------------------------------------------
# Quality control.
is.mito <- grepl("^MT-", rownames(sce))
qcstats <-
    scater::perCellQCMetrics(sce, subsets = list(Mito = is.mito))
filtered <-
    scater::quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
## sacamos genes
sce <- sce[, !filtered$discard]
## single cell experiment

# Normalization.
sce <- scater::logNormCounts(sce)

# Feature selection.
dec <- scran::modelGeneVar(sce)
hvg <- scran::getTopHVGs(dec, prop = 0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- scater::runPCA(sce, ncomponents = 25, subset_row = hvg)
## usamos obj de single cell y solo queremos usar los genes de hvg, por eso el subset row
## es la matrix morada
## esto es para reducir la complejidad

sce <- scater::runUMAP(sce, dimred = 'PCA', external_neighbors = TRUE)
## esto utiliza la matriz morada y es otra manera de reducir dimensiones
## lo puedes usar para visualizar los datos

# Clustering.
g <- scran::buildSNNGraph(sce, use.dimred = 'PCA')
## creamos una red para ver que celulas se parecen a otras
## smallest near neighbor
## usa los datos de pca y crea una matriz

sce$clusters <- factor(igraph::cluster_louvain(g)$membership)


## ----'quick_intro_04'------------------------------------------------------------------------------------------------
# Visualization.
scater::plotUMAP(sce, colour_by = "clusters")
## es un obj de ggplot

## ----'quick_intro_05', eval = FALSE----------------------------------------------------------------------------------
## # Interactive visualization
p <- scater::plotUMAP(sce, colour_by = "clusters")

plotly::ggplotly(p)
## se hace interactivo, pero pesa mucho el doc
