## taken from https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/03-quality-control.R#L57-L320

## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
## Data
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
# Annotate each gene with its chromosome location
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SEQNAME"
)
# Identify the mitochondrial genes
is.mito <- which(location == "MT")

library('scater')
sce.416b <- addPerCellQC(sce.416b,
    subsets = list(Mito = is.mito))
## el resultado de esta funcion es un obj de sce
## estamos crenado nuevas 
## es por cada celula
## la info viene en colData
## detected -> genes con expr mayor a ceropara nuestro obj sce.416b y los alt exp
## percellqc hace el analisis 
#colData()

## ----qc_metrics, cache=TRUE, dependson='all_code'--------------------------------------------------------------------
plotColData(sce.416b, x = "block", y = "detected")
## estas vars deben estar definidad es en colnames(colData(sce.416b))

plotColData(sce.416b, x = "block", y = "detected") +
    scale_y_log10()

plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)


## ----all_code_part2, cache = TRUE, dependson='all_code'--------------------------------------------------------------
# Example thresholds
## definie nuestros thresholds
qc.lib <- sce.416b$sum < 100000 ## lecturos
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10 # pocentaje de expr ercc
qc.mito <- sce.416b$subsets_Mito_percent > 10 ## mas del 10percent en mito
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito
discard <- qc.nexprs | qc.spike | qc.mito
discard <- c.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib),
    NExprs = sum(qc.nexprs),
    SpikeProp = sum(qc.spike),
    MitoProp = sum(qc.mito),
    Total = sum(discard)
)

## cuantas celulas no pasan los thresholds?
## DataFrame with 1 row and 5 columns
##    LibSize    NExprs SpikeProp  MitoProp     Total
##  <integer> <integer> <integer> <integer> <integer>
##   1         3         0        19        14        33
## 14 no pasaron el threshold de mito
## en total deberiamos eliminar 33...

args(isOutlier)
plotColData(sce.416b, x = "block", y = "sum")
plotColData(sce.416b, x = "block", y = "sum") + scale_y_log10()

## scater::isOutliers has a set of assumes

## todos estos son vectores logicos
qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
## en nuestro caso consideramos que solo las bajas no son malas
qc.nexprs2 <- isOutlier(sce.416b$detected, log = TRUE,
    type = "lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extract the thresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
# Summarize the number of cells removed for each reason.
DataFrame(
    LibSize = sum(qc.lib2),
    NExprs = sum(qc.nexprs2),
    SpikeProp = sum(qc.spike2),
    MitoProp = sum(qc.mito2),
    Total = sum(discard2)
)

## More checks
plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)

batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
## mexclamos fenotipo con el bloque
table(batch)
## tenemos 4 

## QC questions

## Was qc.lib necessary for creating discord?
## yes, because in that way we see if a cell was broken 
## yes, porque descarto 2

## > table(qc.spike,qc.mito,qc.lib)
## , , qc.lib = FALSE
##
##       qc.mito
## qc.spike FALSE TRUE
##   FALSE   159   12
##   TRUE     16    2
##
##, , qc.lib = TRUE
##
##        qc.mito
## qc.spike FALSE TRUE
##   FALSE     2    0
##   TRUE      1    0


## Which filter discarded more cells? discard or discard2?
## discard -33, discard2 -6

## By considering the sample batch, did we discard more cells using automatic threshold detection?



qc.lib3 <- isOutlier(sce.416b$sum,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.nexprs3 <- isOutlier(sce.416b$detected,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher",
    batch = batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher",
    batch = batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extract the thresholds
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib3),
    NExprs = sum(qc.nexprs3),
    SpikeProp = sum(qc.spike3),
    MitoProp = sum(qc.mito3),
    Total = sum(discard3)
)


## ----use_case, cache=TRUE, dependson= c('all_code', 'all_code_part2')------------------------------------------------
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)
## addpercel agrega un monton de info


## reminder: no nos gusta tener alto porcentaje de ERCC aka spike-ins
plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")
## las distribuciones de porcentaje de spike-ins que alinearon nos es igual entre donadores
## D10 y D3 presentan muchisimas celulas con alto porcentaje
## D7 todavia pasa de panzaso (-Leo)

hist(sce.grun$altexp_ERCC_percent[sce.grun$donor == 2])
## a ojo, como el 40% de celulas salieron mal
## la distrib no es simetrica, que es algo que habiamos asumido en el modelo

discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor)
## funciona bien pero no para las celulas con distribucion fea d10 y d3

discard.ercc2 <- isOutlier(
    sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor,
    subset = sce.grun$donor %in% c("D17", "D2", "D7")
)
## esta vez nada mas usamos la info de d17, d2 y d2 para poner los limites de lo que aceptas
##      elegimos eso utlimo a ojo, viendo su plot (el de plotColData)
## en esta ocasion si descarto las celulas chafas en d10 y d3
## %in% = 

plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc)
)
plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc2)
)

# Add info about which cells are outliers
sce.416b$discard <- discard2

# Look at this plot for each QC metric
plotColData(
    sce.416b,
    x = "block",
    y = "sum",
    colour_by = "discard",
    other_fields = "phenotype"
) +
    facet_wrap( ~ phenotype) +
    scale_y_log10()

# Another useful diagnostic plot
plotColData(
    sce.416b,
    x = "sum",
    y = "subsets_Mito_percent",
    colour_by = "discard",
    other_fields = c("block", "phenotype")
) +
    facet_grid(block ~ phenotype)


## ----use_case_pbmc, cache=TRUE, dependson='all_code'-----------------------------------------------------------------
library('BiocFileCache')
bfc <- BiocFileCache()
raw.path <-
    bfcrpath(
        bfc,
        file.path(
            "http://cf.10xgenomics.com/samples",
            "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
        )
    )
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library('DropletUtils')
library('Matrix')
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

bcrank <- barcodeRanks(counts(sce.pbmc))
## bar code rank (bcrank)

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(
    bcrank$rank[uniq],
    bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    cex.lab = 1.2
)

## nad mas se estan grafi los ranks unicos
## y-- numero total de UMI Unique Molecular Identifiers
## x -- no de ranksomics
## asi se suelen ver las graficas de 10x gen

## UMI generalmnete es como el no unico de genes

abline(h = metadata(bcrank)$inflection,
    col = "darkgreen",
    lty = 2)
abline(h = metadata(bcrank)$knee,
    col = "dodgerblue",
    lty = 2)
legend(
    "bottomleft",
    legend = c("Inflection", "Knee"),
    col = c("darkgreen", "dodgerblue"),
    lty = 2,
    cex = 1.2
) 
## el metodo estadistico va a querer determinar entre el intervalo de las dos lineas
##      que droplets tienen celulas y cuales no

set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
## empty Drops tinee el metodo estadistico 
##          que se utiliza para decidir si un droplet esta vacio o no

# See ?emptyDrops for an explanation of why there are NA # values.
summary(e.out$FDR <= 0.001)

set.seed(100)
limit <- 100
all.out <-
    emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE)
# Ideally, this histogram should look close to uniform.
# Large peaks near zero indicate that barcodes with total
# counts below 'lower' are not ambient in origin.
hist(all.out$PValue[all.out$Total <= limit &
        all.out$Total > 0],
    xlab = "P-value",
    main = "",
    col = "grey80")

sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
sce.pmbc <- addPerCellQC(sce.pbmc, subsets = list(MT = is.mito))
discard.mito <-
    isOutlier(sce.pmbc$subsets_MT_percent, type = "higher")
plot(
    sce.pmbc$sum,
    sce.pmbc$subsets_MT_percent,
    log = "x",
    xlab = "Total count",
    ylab = "Mitochondrial %"
)
abline(h = attr(discard.mito, "thresholds")["higher"], col = "red")


## ----marking, cache=TRUE, dependson='use_case'-----------------------------------------------------------------------
# Removing low-quality cells
# Keeping the columns we DON'T want to discard
filtered <- sce.416b[,!discard2]
# Marking low-quality cells
marked <- sce.416b
marked$discard <- discard2


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()
