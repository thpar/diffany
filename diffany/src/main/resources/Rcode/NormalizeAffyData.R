expressionSet <- rma(rawProbeData);
probesets <- featureNames(expressionSet);
samples <- sampleNames(expressionSet);
expressionMatrix <- exprs(expressionSet);
