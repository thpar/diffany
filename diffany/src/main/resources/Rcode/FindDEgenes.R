#source("http://bioconductor.org/biocLite.R");
#biocLite("limma");
#install.packages("statmod") 

library(limma);
library(statmod);

# GENERIC COMPARISON BETWEEN ALL STRESS AND ALL CONTROL CONDITIONS (COARSE-GRAINED)
design_stress <- model.matrix(~0+Stress, data=pData(expressionSet));
colnames(design_stress) <- gsub("Stress","",colnames(design_stress));

# CHOSE WHETHER OR NOT TO MODEL THE REPLICATES AS RANDOM VARIABLES
corfit_stress <- duplicateCorrelation(expressionSet, design=design_stress, block=targets$Replicate);
fit_stress <- lmFit(expressionSet, design_stress, block=targets$Replicate, cor=corfit_stress$consensus);
#fit_stress <- lmFit(expressionSet, design_stress);

contrasts_stress <- makeContrasts(StressvsControl=mannitol-control, levels=design_stress);
fit2_stress <- contrasts.fit(fit_stress, contrasts_stress);
efit_stress <- eBayes(fit2_stress);
toptable_stress <- topTable(efit_stress, coef="StressvsControl", number=Inf);
topIDs_stress <- row.names(toptable_stress);

# SPECIFIC COMPARISONS BETWEEN TIME-DEPENDENT STRESS AND CONTROL EXPERIMENTS
design_stress_time <- model.matrix(~0+Setup, data=pData(expressionSet));
colnames(design_stress_time) <- gsub("Setup","",colnames(design_stress_time));

# CHOSE WHETHER OR NOT TO MODEL THE REPLICATES AS RANDOM VARIABLES
corfit_stress_time <- duplicateCorrelation(expressionSet, design=design_stress_time, block=targets$Replicate);
fit_stress_time <- lmFit(expressionSet, design_stress_time, block=targets$Replicate, cor=corfit_stress_time$consensus);
#fit_stress_time <- lmFit(expressionSet, design_stress_time);

# TURN THIS ON WHEN ALL CONTROL MEASURES ARE GROUPED
# contrasts_stress_time <- makeContrasts(response_3h=mannitol_3h-control_0h, response_12h=mannitol_12h-control_0h, response_24h=mannitol_24h-control_0h, response_1.5h=mannitol_1.5h-control_0h, levels=design_stress_time);

contrasts_stress_time <- makeContrasts(response_3h=mannitol_3h-control_3h, response_12h=mannitol_12h-control_12h, response_24h=mannitol_24h-control_24h, response_1.5h=mannitol_1.5h-control_1.5h, levels=design_stress_time);

fit2_stress_time <- contrasts.fit(fit_stress_time, contrasts_stress_time);
efit_stress_time <- eBayes(fit2_stress_time);

#toptable_stress_time_3 <- topTable(efit_stress_time, coef="response_3h", number=Inf, p.value=0.05);
toptable_stress_time_3 <- topTable(efit_stress_time, coef="response_3h", number=Inf);
topIDs_stress_time_3 <- row.names(toptable_stress_time_3);

toptable_stress_time_12 <- topTable(efit_stress_time, coef="response_12h", number=Inf);
topIDs_stress_time_12 <- row.names(toptable_stress_time_12);

toptable_stress_time_24 <- topTable(efit_stress_time, coef="response_24h", number=Inf);
topIDs_stress_time_24 <- row.names(toptable_stress_time_24);

toptable_stress_time_1.5 <- topTable(efit_stress_time, coef="response_1.5h", number=Inf);
topIDs_stress_time_1.5 <- row.names(toptable_stress_time_1.5);

