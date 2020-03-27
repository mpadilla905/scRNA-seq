## copied source code from: https://github.com/lcolladotor/osca_LIIGH_UNAM_2020/blob/master/04-normalization.R#L57-L140

## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
library('scRNAseq')
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)


# Quality control
library('scater')
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
## cuales genes son mitocondrial
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
    quickPerCellQC(stats,
        percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))
## QUCKKper cell corre de fondo isOutlier

sce.zeisel <- sce.zeisel[, !qc$discard]
## descartamos las celulas de baja calidad 190 de 3000
table(qc$discard)



## ----all_code2, cache=TRUE, dependson='all_code'---------------------------------------------------------------------
# Library size factors
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
## esta funcion nos va a dar el size factor dependiendo de la libreria
## librery size = sum of counts across all genes for each cell

class(lib.sf.zeisel)
length(lib.sf.zeisel)
mean(lib.sf.zeisel)
head(lib.sf.zeisel)

# Examine distribution of size factors
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab = "Log10[Size factor]", col = "grey80")

head(counts(sce.zeisel))

## contar manualmente el library size
ls.zeisel <- colSums(counts(sce.zeisel))
head(ls.zeisel)

identical(lib.sf.zeisel,ls.zeisel)

plot(
    ls.zeisel,
    lib.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)

##########################################################################
## EXERCISE 1

## Are ls.zeisel and lib.sf.zeisel identical?
m = match(names(lib.sf.zeisel),names(ls.zeisel))
ls.zeisel = ls.zeisel[m]
identical(lib.sf.zeisel,ls.zeisel)
## FALSE

## Are they proportional? yes

length(lib.sf.zeisel) ## 2815
length(ls.zeisel) ## 2815

sum(lib.sf.zeisel/ls.zeisel)/length(lib.sf.zeisel) # 6.686916e-05
x = lib.sf.zeisel/ls.zeisel
## checar si en x todo es igual
## plot(lib.sf.zeisel, ls.zeisel) ## ninguno se sale de la diagonal perf

# otra opcion
table(lib.sf.zeisel/ls.zeisel)
## 6.68691586731496e-05 
##                2815 

## Compute lib.sf.zeisel manually

##       Library size = sum of counts across all genes for each cell
##       The "library size factor" for each cell is then scaled so that 
##             the average size factor across all cells is 1

## First compute the sums
zeisel_sums <- colSums(counts(sce.zeisel))
identical(zeisel_sums, ls.zeisel)
## TRUE

## Next, make them have unity mean
zeisel_size_factors <- zeisel_sums/mean(zeisel_sums)
identical(zeisel_size_factors, lib.sf.zeisel)
## TRUE

##########################################################################

## ----exercise_solution, cache=TRUE, dependson='all_code'-------------------------------------------------------------
## First compute the sums
zeisel_sums <- colSums(counts(sce.zeisel))
identical(zeisel_sums, ls.zeisel)

## Next, make them have unity mean
zeisel_size_factors <- zeisel_sums/mean(zeisel_sums)
identical(zeisel_size_factors, lib.sf.zeisel)


## ----all_code3, cache=TRUE, dependson='all_code2'--------------------------------------------------------------------
# Normalization by convolution

library('scran')
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)
# Compute deconvolution size factors
deconv.sf.zeisel <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

# Examine distribution of size factors
summary(deconv.sf.zeisel)
hist(log10(deconv.sf.zeisel), xlab = "Log10[Size factor]",
    col = "grey80")
plot(
    ls.zeisel,
    deconv.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)


## ----all_code4, cache=TRUE, dependson='all_code3'--------------------------------------------------------------------
# Library size factors vs. convolution size factors

# Colouring points using the supplied cell-types
plot(
    lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = 'xy',
    pch = 16,
    col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()
