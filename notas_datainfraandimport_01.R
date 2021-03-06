## copied from https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/02-data-infrastructure-and-import.R

## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
## datos de un exp


# Load the SingleCellExperiment package
library('SingleCellExperiment')


# Extract the count matrix from the 416b dataset
counts.416b <- counts(sce.416b)
# Construct a new SCE from the counts matrix
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Inspect the object we just created
sce
# class: SingleCellExperiment
# dim: 46604 192
# metadata(2): favourite_genes analyst
# assays(2): counts logcounts
# rownames(46604): ENSMUSG00000102693 ENSMUSG00000064842 ... ENSMUSG00000095742 CBFB-MYH11-mcherry
# rowData names(4): Length mean detected chromosome
# colnames(192): SLX-9555.N701_S502.C89V9ANXX.s_1.r_1 SLX-9555.N701_S503.C89V9ANXX.s_1.r_1 ...
# SLX-11312.N712_S508.H5H5YBBXX.s_8.r_1 SLX-11312.N712_S517.H5H5YBBXX.s_8.r_1
# colData names(23): Source Name cell line ... total sizeFactor
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(2): ERCC SIRV

## How big is it? 40 MB
pryr::object_size(sce)

dim(counts.416b)
# [1] 46604   192

# Access the counts matrix from the assays slot
# WARNING: This will flood RStudio with output!

# 1. The general method
assay(sce, "counts")[1:6, 1:3]
## nos va crear la matriz

# 2. The special method for the assay named "counts"
counts(sce)[1:6, 1:3]


## vaos a normalizar las cuentas y crea la tabla correspondiente, otra tabla azul
sce <- scater::logNormCounts(sce)
# Inspect the object we just updated
sce

## How big is it? ahora es de 112MB
pryr::object_size(sce) # ya estaba y agregando la tabla nueva
## la tabla es la que dice logcounts en assays

## sce se reescribio pero tomando la info que

# 1. The general method
assay(sce, "logcounts")[1:6, 1:3]
# 2. The special method for the assay named "logcounts"
logcounts(sce)[1:6, 1:3]

# assign a new entry to assays slot
assay(sce, "counts_100") <- assay(sce, "counts") + 100
# List the assays in the object
assays(sce) ### diferentes matrices en nuestro obj
assayNames(sce) # esto nada mas nos ensena el nombre de nuestras matrices de cuentas

## How big is it?
pryr::object_size(sce)

# ya tenemos las cuentas, necesitamos ahora la info de las celulas, del experiment
## para eso tenemos la tabla de Cell Metada y ahora cada col la llamamos con colData

# Extract the sample metadata from the 416b dataset
colData.416b <- colData(sce.416b)
    ## colData es data.frame de Bioconductor pero mas flexible
    ## aqui tenemos la info de fenotipo

# Add some of the sample metadata to our SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Inspect the object we just updated
sce
# Access the sample metadata from our SCE
colData(sce)
# Access a specific column of sample metadata from our SCE
table(sce$block)

# Example of function that adds extra fields to colData
sce <- scater::addPerCellQC(sce.416b)
## tenemos la info de antes mas otras 12 cols

# Access the sample metadata from our updated SCE
colData(sce)

# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

## Add the lognorm counts again
sce <- scater::logNormCounts(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just wild type cells
# Remember, cells are columns of the SCE
sce[, sce$phenotype == "wild type phenotype"]
x = sce$phenotype == "wild type phenotype"
## x es un vector logico tam ncol(sce)


# Access the feature metadata from our SCE
# It's currently empty!
rowData(sce)
## esto es la primera tablita, de metadata del gen

# Example of function that adds extra fields to rowData
sce <- scater::addPerFeatureQC(sce)
## nos va a adar info sobre cada gen

# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)


# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
    # Annotation hub tiene un monton de datos genomicos
query(ah, c("Mus musculus", "Ensembl", "v97"))
    # pedimos especificamente esta info del anootation hub

# Annotate each gene with its chromosome location
ensdb <- ah[["AH73905"]]
    ## el doble corchete es igual a
    ## [["id"]] descargar datos y cargarlos en la sesion de R
## download_data = function(ah, id) { access_data( sh, id) }

chromosome <- mapIds(ensdb,
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
 ## de ensdb tomamos ...

rowData(sce)$chromosome <- chromosome

# Access the feature metadata from our updated SCE
rowData(sce)

## How big is it?
pryr::object_size(sce)

# E.g., subset data to just genes on chromosome 3
# NOTE: which() needed to cope with NA chromosome names
sce[which(rowData(sce)$chromosome == "3"), ]

# Access the metadata from our SCE
# It's currently empty!
metadata(sce)

# The metadata slot is Vegas - anything goes
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
    analyst = c("Pete"))

# Access the metadata from our updated SCE
metadata(sce)

# E.g., add the PCA of logcounts
# NOTE: We'll learn more about PCA later
sce <- scater::runPCA(sce)
# Inspect the object we just updated
sce
# Access the PCA matrix from the reducedDims slot
reducedDim(sce, "PCA")[1:6, 1:3]

# E.g., add a t-SNE representation of logcounts
# NOTE: We'll learn more about t-SNE later
sce <- scater::runTSNE(sce)
# Inspect the object we just updated
sce
# Access the t-SNE matrix from the reducedDims slot
head(reducedDim(sce, "TSNE"))

# E.g., add a 'manual' UMAP representation of logcounts
# NOTE: We'll learn more about UMAP later and a
# 		  simpler way to compute it.
u <- uwot::umap(t(logcounts(sce)), n_components = 2)
# Add the UMAP matrix to the reducedDims slot
# Access the UMAP matrix from the reducedDims slot
reducedDim(sce, "UMAP") <- u

# List the dimensionality reduction results stored in # the object
reducedDims(sce)

# Extract the ERCC SCE from the 416b dataset
ercc.sce.416b <- altExp(sce.416b, "ERCC")
    # altExp: info de experimentos alternaticos, e.g.
    # ercc es un tipo de spike-ins, son 92 seqs para evaluar
    # preparacion de pcr de diferentes longs
# Inspect the ERCC SCE
ercc.sce.416b

# Add the ERCC SCE as an alternative experiment to our SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Inspect the object we just updated
sce

## How big is it?
pryr::object_size(sce)

# List the alternative experiments stored in the object
altExps(sce)

# Subsetting the SCE by sample also subsets the
# alternative experiments
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## How big is it?
pryr::object_size(sce.subset)

# Extract existing size factors (these were added
# when we ran scater::logNormCounts(sce))
head(sizeFactors(sce))
## agrgaos info del numero de genes

# 'Automatically' replace size factors
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# 'Manually' replace size factors
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))

## ----ercc_exercise, cache = TRUE, dependson='all_code'---------------------------------------------------------------
## Leer tabla https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt con read.delim()
## Read the data from the web
ercc_info <-
    read.delim(
        'https://tools.thermofisher.com/content/sfs/manuals/cms_095046.txt',
        as.is = TRUE,
        row.names = 2,
        check.names = FALSE
    )

## Explorando el objecto de altExp(sce, 'ERCC')
dim(altExp(sce, 'ERCC')) # 92 genes x celulas
dim(ercc_info) ## info para cada gene (92)


assayNames(altExp(sce, 'ERCC'))
dim(assay(altExp(sce, 'ERCC'), 'counts'))


## Usar los ERCC ID para alinear esta tabla con el objeto sce (ERCC alt experiment)

## Match the ERCC data
m <- match(rownames(altExp(sce, "ERCC")), rownames(ercc_info))
ercc_info <- ercc_info[m, ]

## Normalize the ERCC counts
#altExp(sce, "ERCC") <- scater::logNormCounts(altExp(sce, "ERCC"))
## mejor no xd


## Usar plot() para graficar concentration in Mix 1 (attomoles/ul) vs las cuentas de ERCC de nuestro sce(en alt exp)
i <- 1
assay(altExp(sce, 'ERCC'), 'counts')[, i]
plot(ercc_info[,"concentration in Mix 1 (attomoles/ul)"] ~ counts(altExp(sce, "ERCC"))[,i] )



## ----ercc_solution_plots, cache = TRUE, dependson='ercc_exercise'----------------------------------------------------
for (i in seq_len(2)) {
    plot(
        log2(10 * ercc_info[, "concentration in Mix 1 (attomoles/ul)"] + 1) ~
            log2(counts(altExp(sce, "ERCC"))[, i] +
                    1),
        xlab = "log norm counts",
        ylab = "Mix 1: log2(10 * Concentration + 1)",
        main = colnames(altExp(sce, "ERCC"))[i],
        xlim = c(min(logcounts(
            altExp(sce, "ERCC")
        )), max(logcounts(
            altExp(sce, "ERCC")
        )))
    )
    abline(0, 1, lty = 2, col = 'red')
}


## de hacer un pdf el chiste es que puedes explorar todas las celulas y ver si hay un problema en especifico
## ver patrones de errores
pdf('ERCC_example.pdf')
for (i in seq_len(ncol(sce))) {
    message(paste(Sys.time(), 'plotting cell', i))
    plot(
        log2(10 * ercc_info[, "concentration in Mix 1 (attomoles/ul)"] + 1) ~
            log2(counts(altExp(sce, "ERCC"))[, i] +
                    1),
        xlab = "log norm counts",
        ylab = "Mix 1: log2(10 * Concentration + 1)",
        main = colnames(altExp(sce, "ERCC"))[i],
        xlim = c(min(logcounts(
            altExp(sce, "ERCC")
        )), max(logcounts(
            altExp(sce, "ERCC")
        )))
    )
    abline(0, 1, lty = 2, col = 'red')
}
dev.off()
## en el pdf vamos a tener, una grafica por celula,
## de la concentrancion de cada ERCC gene contra las lecturas de cada ERCC gen


## ----all_code_part2, cache=TRUE--------------------------------------------------------------------------------------
# Download example data processed with CellRanger
# Aside: Using BiocFileCache means we only download the
#        data once
library('BiocFileCache')
bfc <- BiocFileCache()
pbmc.url <-
    paste0(
        "http://cf.10xgenomics.com/samples/cell-vdj/",
        "3.1.0/vdj_v1_hs_pbmc3/",
        "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
    )
pbmc.data <- bfcrpath(bfc, pbmc.url)

# Extract the files to a temporary location
untar(pbmc.data, exdir = tempdir())

# List the files we downloaded and extracted
# These files are typically CellRanger outputs
pbmc.dir <- file.path(tempdir(),
    "filtered_feature_bc_matrix")
list.files(pbmc.dir)

# Import the data as a SingleCellExperiment
library('DropletUtils')
sce.pbmc <- read10xCounts(pbmc.dir)
# Inspect the object we just constructed
sce.pbmc

## How big is it?
pryr::object_size(sce.pbmc)

# Store the CITE-seq data in an alternative experiment
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Inspect the object we just updated
sce.pbmc

## How big is it?
pryr::object_size(sce.pbmc)

# Download example data processed with scPipe
library('BiocFileCache')
bfc <- BiocFileCache()
sis_seq.url <-
    "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extract the files to a temporary location
unzip(sis_seq.data, exdir = tempdir())

# List (some of) the files we downloaded and extracted
# These files are typical scPipe outputs
sis_seq.dir <- file.path(tempdir(),
    "SIS-seq_script-master",
    "data",
    "BcorKO_scRNAseq",
    "RPI10")
list.files(sis_seq.dir)

# Import the data as a SingleCellExperiment
library('scPipe')
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Inspect the object we just constructed
sce.sis_seq

## How big is it?
pryr::object_size(sce.sis_seq)

# Download example bunch o' files dataset
library('BiocFileCache')
bfc <- BiocFileCache()
lun_counts.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
    )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <-
    paste0("https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt")
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extract the counts files to a temporary location
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

# List the files we downloaded and extracted
list.files(lun_counts.dir)

# Import the count matrix (for 1 plate)
lun.counts <- read.delim(
    file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
# Store the gene lengths for later
gene.lengths <- lun.counts$Length
# Convert the gene counts to a matrix
lun.counts <- as.matrix(lun.counts[, -1])

# Import the sample metadata
lun.coldata <- read.delim(lun_coldata.data,
    check.names = FALSE,
    stringsAsFactors = FALSE)
library('S4Vectors')
lun.coldata <- as(lun.coldata, "DataFrame")

# Match up the sample metadata to the counts matrix
m <- match(colnames(lun.counts),
    lun.coldata$`Source Name`)
lun.coldata <- lun.coldata[m,]

# Construct the feature metadata
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construct the SingleCellExperiment
lun.sce <- SingleCellExperiment(
    assays = list(assays = lun.counts),
    colData = lun.coldata,
    rowData = lun.rowdata
)
# Inspect the object we just constructed
lun.sce

## How big is it?
pryr::object_size(lun.sce)


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()

## SCE OBJECT EXERCISES

## which function defines the sce class?
## SingleCellExperiment::SingleCellExperiment

## what are the minimum type of tables an sce object contains?
## infoGenes : rowData()
## number of reads overlapping each gene for each cell: assays
## info about cells : colData()
## optionally: PCA, TSNE, alternative experiments
## random info (metadata)

## where are the colnames(sce) used?
head(colnames(sce))
## column names of the assays + rownames of the colData, nombres de las celulas

## Similarly, where are the rownames(sce) used?
## head(rownames(sce))
## head(rownames(rowData(sce)))
## nombres de los renglones de la matriz azul, o de la matriz verde


## How many principal components did we compute?
## dim(reducedDim(sce, "PCA"))
## dim(sce)
## head(reducedDim(sce, "PCA")) 50 PCs

## Which three chromosomes have the highest mean gene expression?
rowData(sce)
sort(with(rowData(sce), tapply(mean, chromosome, base::mean)), decreasing = TRUE )





