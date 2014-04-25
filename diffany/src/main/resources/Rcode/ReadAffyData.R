#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");

library(affy);

pheno <- read.AnnotatedDataFrame("RNAtargets_Sofie.txt");
rawProbeData <- ReadAffy(phenoData=pheno);
#probeValues <- exprs(rawProbeData);
#probeNames <- probeNames(rawProbeData);
