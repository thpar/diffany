#source("http://bioconductor.org/biocLite.R");
#biocLite("limma");

library(limma);

design_stress <- model.matrix(~0+Stress, data=pData(expressionSet));
fit_stress <- lmFit(expressionSet, design_stress);
contrasts_stress <- makeContrasts(StressvsControl=Stressmannitol25-Stresscontrol, levels=design_stress);
fit2_stress <- contrasts.fit(fit_stress, contrasts_stress);
efit_stress <- eBayes(fit2_stress);
toptable_stress <- topTable(efit_stress, coef=1, number=Inf, p.value=0.05);
topIDs_stress <- row.names(toptable_stress);

#design_stress_time <- model.matrix(~Time*Stress, data=pData(expressionSet));
#colnames(design_stress_time) <- gsub(":","_",colnames(design_stress_time));
#colnames(design_stress_time)[1] <- Intercept;

#contrasts_stress_time <- makeContrasts(Time3h_Stressmannitol25-Time3h, Time12h_Stressmannitol25-Time12h, Time24h_Stressmannitol25-Time24h, levels=design_stress_time);
#fit_stress_time <- lmFit(expressionSet, design_stress_time);
#fit2_stress_time <- contrasts.fit(fit_stress_time, contrasts_stress_time);
#efit_stress_time <- eBayes(fit2_stress_time);

#toptable_stress_time_3 <- topTable(efit_stress_time, coef=1, number=Inf);
#topIDs_stress_time_3 <- row.names(toptable_stress_time_3);

#toptable_stress_time_12 <- topTable(efit_stress_time, coef=2, number=Inf);
#topIDs_stress_time_12 <- row.names(toptable_stress_time_12);

#toptable_stress_time_24 <- topTable(efit_stress_time, coef=3, number=Inf);
#topIDs_stress_time_24 <- row.names(toptable_stress_time_24);

