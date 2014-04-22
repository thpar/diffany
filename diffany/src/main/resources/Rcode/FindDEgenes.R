#source("http://bioconductor.org/biocLite.R");
#biocLite("limma");

library(limma);

#combn <- factor(paste(pData(pheno)[, 1], pData(pheno)[, 2], sep="_"));
#design <- model.matrix(~combn);

#design <- cbind(c=1,mutvsc=c(0,0,0,1,1,1)); 

conditions <- factor(c("control","control","control","mannitol","mannitol","mannitol","control","control","control","mannitol","mannitol","mannitol","control","control","control","mannitol","mannitol","mannitol","control","control","control","mannitol","mannitol","mannitol"));
times <- factor(c(3,3,3,3,3,3,12,12,12,12,12,12,24,24,24,24,24,24,1.5,1.5,1.5,1.5,1.5,1.5));
design <- model.matrix(~times*conditions);

fit <- lmFit(expressionSet, design);
efit <- eBayes(fit);
toptable <- topTable(efit, coef=2, number=50);
topIDs <- row.names(toptable);
