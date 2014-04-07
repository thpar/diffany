#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");
#biocLite("IRanges");

library(affy);
library(Biobase);

pheno <- read.AnnotatedDataFrame("RNAtargets_Sofie.txt");
rawProbeData <- ReadAffy(phenoData=pheno);
#probeValues <- exprs(rawProbeData);
#probeNames <- probeNames(rawProbeData);
