#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("limma");

library(affy);
library(limma);

targets <- readTargets("RNAtargets_Sofie.txt");

# WHEN TURNING THIS ON, ALL CONTROL MEASURES ARE GROUPED
#targets[targets$Stress=="control",]$Time <- "0h";

cols <- c("Stress","Time");
targets[, "Setup"] <- do.call(paste, c(targets [cols], sep="_"));
pheno <- new("AnnotatedDataFrame", data=targets);

rawProbeData <- ReadAffy(phenoData=pheno);
#probeValues <- exprs(rawProbeData);
#probeNames <- probeNames(rawProbeData);
