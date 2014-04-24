#source("http://bioconductor.org/biocLite.R");
#biocLite("limma");

library(limma);

design_stress <- model.matrix(~Stress, data=pData(expressionSet))
design_stress_time <- model.matrix(~Time*Stress, data=pData(expressionSet))

fit_stress <- lmFit(expressionSet, design_stress);
fit_stress_time <- lmFit(expressionSet, design_stress_time);

efit_stress <- eBayes(fit_stress);
efit_stress_time <- eBayes(fit_stress_time);

toptable_stress <- topTable(efit_stress, coef=2, number=Inf);
toptable_stress_time <- topTable(efit_stress_time, coef=2, number=Inf);

topIDs_stress <- row.names(toptable_stress);
topIDs_stress_time <- row.names(toptable_stress_time);
