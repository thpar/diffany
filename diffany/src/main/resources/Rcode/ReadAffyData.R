#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("IRanges");
#biocLite("limma");
#biocLite("affyPLM");
#biocLite("org.Dm.eg.db");
library(affy);
library(limma);
library(affyPLM);
library(org.Dm.eg.db);
library(Biobase)
pheno <- read.AnnotatedDataFrame("RNAtargets.txt");
rawProbeData <- ReadAffy(phenoData=pheno);
#probeValues <- exprs(rawProbeData);
#probeNames <- probeNames(rawProbeData);
