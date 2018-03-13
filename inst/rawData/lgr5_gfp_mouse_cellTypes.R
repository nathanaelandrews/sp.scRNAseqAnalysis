#run from the package root as: source('inst/rawData/lgr5_gfp_mouse_cellTypes.R')

library(sp.scRNAseq)
library(tidyverse)

#ideally you would only run this if the operating system is windows
#there is some info here: https://stackoverflow.com/questions/4747715/how-to-check-the-os-within-r
#but I can't check it without a windows computer
#You should be able to write something like:
if(Sys.info()['sysname'] == "Windows") {extrafont::loadfonts(device="win")}
#extrafont::loadfonts(device="win")

#load("counts_raw_and_log_180112.rda")
#load("counts_ercc.rda")
#this is better because it will work for anyone that has your package installed
#and therefore is more reproduciable.
load(system.file('rawData/counts_raw_and_log_180112.rda', package = "sp.scRNAseqAnalysis"))
load(system.file('rawData/counts_ercc.rda', package = "sp.scRNAseqAnalysis"))

#remove .htseq suffix
colnames(counts) <- str_replace(colnames(counts), "(.*)\\.htseq$", "\\1")
colnames(counts.ercc) <- str_replace(colnames(counts.ercc), "(.*)\\.htseq$", "\\1")

#Select Singlets
s <- grepl("SRR", colnames(counts)) | grepl("Singlet", colnames(counts)) | grepl("NJA00102", colnames(counts)) | grepl("NJA00103", colnames(counts)) | grepl("NJA00104", colnames(counts)) | grepl("NJA00109", colnames(counts)) | grepl("NJA00204", colnames(counts)) | grepl("NJA00205", colnames(counts)) | grepl("NJA00206", colnames(counts))

counts.ercc[, grepl("SRR", colnames(counts.ercc))] <- NA
cObjSng <- spCounts(counts[, s], counts.ercc[, s])
cObjMul <- spCounts(counts[, !s], counts.ercc[, !s])

uObj <- spUnsupervised(cObjSng, seed = 877239)

classification <- getData(uObj, "classification")
names(classification) <- getData(uObj, "tsne") %>%
  rownames()

SI.Stem <- names(classification)[classification %in% c("A1", "J1", "K1", "N1", "T1")]
SI.TA.Cells <- names(classification)[classification %in% c("R1", "C1", "B1")]
SI.Enterocyte.Progenitor <- names(classification)[classification %in% c("E1", "Q1")]
SI.Early.Enterocyte <- names(classification)[classification == "G1"]
SI.Late.Enterocyte <- names(classification)[classification == "M1"]
SI.Paneth <- names(classification)[classification %in% c("P1", "O1")]
SI.Goblet <- names(classification)[classification %in% c("F1", "L1")]
Blood.Cells <- names(classification)[classification == "D1"]
SI.Tuft.Cells <- names(classification)[classification == "H1"]
SI.Enteroendocrine <- names(classification)[classification == "S1"]
Colon <- names(classification)[classification == "I1"]

#tSNE for colon cells
cs <- colnames(counts) %in% Colon

#Colon object
colcObjSng <- spCounts(counts[, cs], counts.ercc[, cs])

#tSNE for Colon cell and expression of markers
C.uObj <- spUnsupervised(colcObjSng)

#Select and rename clusters according to cell types in colon
C.classification <- getData(C.uObj, "classification")
names(C.classification) <- getData(C.uObj, "tsne") %>%
  rownames()

C.Stem <- names(C.classification)[C.classification == "B1"]
C.Goblet <- names(C.classification)[C.classification == "C1"]
C.Colonocytes <- names(C.classification)[C.classification == "A1"]

#Organize new classifications and insert them into uObj which contains all cells
celltypes <- data.frame(
  Samples = c(
    SI.Early.Enterocyte, SI.Enterocyte.Progenitor, SI.Enteroendocrine,
    SI.Goblet, SI.Late.Enterocyte, SI.Paneth, SI.Stem, SI.TA.Cells,
    C.Colonocytes, C.Goblet, C.Stem, SI.Tuft.Cells, Blood.Cells
  ),
  class = c(
    rep("SI.Early.Enterocyte", length(SI.Early.Enterocyte)),
    rep("SI.Enterocyte.Progenitor", length(SI.Enterocyte.Progenitor)),
    rep("SI.Enteroendocrine", length(SI.Enteroendocrine)),
    rep("SI.Goblet", length(SI.Goblet)),
    rep("SI.Late.Enterocyte", length(SI.Late.Enterocyte)),
    rep("SI.Paneth", length(SI.Paneth)),
    rep("SI.Stem", length(SI.Stem)),
    rep("SI.TA.Cells", length(SI.TA.Cells)),
    rep("C.Colonocytes", length(C.Colonocytes)),
    rep("C.Goblet", length(C.Goblet)),
    rep("C.Stem", length(C.Stem)),
    rep("SI.Tuft.Cells", length(SI.Tuft.Cells)),
    rep("Blood.Cells", length(Blood.Cells))
  )
)

newclass <- celltypes[match(rownames(getData(uObj, "tsne")), celltypes$Samples),]
newclass <- as.vector(newclass$class)

lgr5_gfp_mouse_cellTypes.newCelltypes <- celltypes

uObj@classification <- newclass
tsneMeans(uObj) <- tsneGroupMeans(getData(uObj, "tsne"), getData(uObj, "classification"))
groupMeans(uObj) <- averageGroupExpression(cObjSng, getData(uObj, "classification"), weighted=FALSE)

#Set Uncertainty to 0 and attempt to identify cell composition of multiplets

uncertainty(uObj) <- rep(0, length(getData(uObj, "uncertainty")))
lgr5_gfp_mouse_cellTypes.uObj <- uObj
save(lgr5_gfp_mouse_cellTypes.newCelltypes, lgr5_gfp_mouse_cellTypes.uObj, file = "data/lgr5_gfp_mouse_cellTypes.rda", compress = "bzip2")



