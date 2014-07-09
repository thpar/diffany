#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("limma");

library(affy);
library(limma);

targets <- readTargets("RNAtargets_Sofie.txt");

cols <- c("Stress","Time");
targets[, "Setup"] <- do.call(paste, c(targets [cols], sep="_"));
pheno <- new("AnnotatedDataFrame", data=targets);

rawProbeData <- ReadAffy(phenoData=pheno, filenames=pheno$FileName);
