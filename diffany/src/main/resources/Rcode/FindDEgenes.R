#source("http://bioconductor.org/biocLite.R");
#biocLite("limma");
#install.packages("statmod"); 

library(limma);
library(statmod);

# SPECIFIC COMPARISONS BETWEEN TIME-DEPENDENT STRESS AND CONTROL EXPERIMENTS
design_stress_time <- model.matrix(~0+Setup, data=pData(expressionSet));
colnames(design_stress_time) <- gsub("Setup","",colnames(design_stress_time));

# CHOSE WHETHER OR NOT TO MODEL THE REPLICATES AS RANDOM VARIABLES
corfit_stress_time <- duplicateCorrelation(expressionSet, design=design_stress_time, block=targets$Replicate);
fit_stress_time <- lmFit(expressionSet, design_stress_time, block=targets$Replicate, cor=corfit_stress_time$consensus);
#fit_stress_time <- lmFit(expressionSet, design_stress_time);

contrasts_stress_time <- makeContrasts(response_3h=mannitol_3h-control_3h, response_12h=mannitol_12h-control_12h, response_24h=mannitol_24h-control_24h, response_1.5h=mannitol_1.5h-control_1.5h, levels=design_stress_time);

fit2_stress_time <- contrasts.fit(fit_stress_time, contrasts_stress_time);
efit_stress_time <- eBayes(fit2_stress_time);

#toptable_3 <- topTable(efit_stress_time, coef="response_3h", number=Inf, p.value=0.05);
toptable_3 <- topTable(efit_stress_time, coef="response_3h", number=Inf);
topIDs_3 <- row.names(toptable_3);

toptable_12 <- topTable(efit_stress_time, coef="response_12h", number=Inf);
topIDs_12 <- row.names(toptable_12);

toptable_24 <- topTable(efit_stress_time, coef="response_24h", number=Inf);
topIDs_24 <- row.names(toptable_24);

toptable_1.5 <- topTable(efit_stress_time, coef="response_1.5h", number=Inf);
topIDs_1.5 <- row.names(toptable_1.5);

