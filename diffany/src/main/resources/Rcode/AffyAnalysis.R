library(affy);
setwd("D:/diffany-osmotic/data-marieke/short-term-osmotic-stress");
ourdata <- ReadAffy();
ed <- exprs(ourdata);
probes <- featureNames(ourdata);