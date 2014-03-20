source("http://bioconductor.org/biocLite.R");
biocLite("affy");
biocLite("affyPLM");
biocLite("org.Dm.eg.db");
library(affy);
library(affyPLM);
library(org.Dm.eg.db);
setwd("D:/diffany-osmotic/data-marieke/short-term-osmotic-stress");
ourdata <- ReadAffy();
ed <- exprs(ourdata);
probes <- featureNames(ourdata);
samp <- sampleNames(ourdata)