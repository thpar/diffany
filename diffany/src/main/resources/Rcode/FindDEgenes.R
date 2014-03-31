#combn <- factor(paste(pData(pheno)[, 1], pData(pheno)[, 2], sep = "_"));
#design <- model.matrix(~combn);
design <- cbind(c=1,mutvsc=c(0,0,0,1,1,1));
fit <- lmFit(expressionSet, design);
efit <- eBayes(fit);
toptable <- topTable(efit, coef = 2);
topIDs <- row.names(toptable);
