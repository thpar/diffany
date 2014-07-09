#source("http://bioconductor.org/biocLite.R");
#biocLite("affy");

library(affy);

expressionSet <- rma(rawProbeData);
probesets <- featureNames(expressionSet);
samples <- sampleNames(expressionSet);
