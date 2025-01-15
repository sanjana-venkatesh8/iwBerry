% test script for Foldiak model
foldiakParams = StimuliParams(nX=8, nY=8, nOrient=4, barLength=8, scanLength=8, ...
    nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=400, ...
    nDendRecord=3, recordingTimeInterval=1, isFoldiak=true);
% stimuliParamsObj = StimuliParams(8, 8, 4, 8, 8, 0, 0, 0, 0, 1000, 400, 3, 1, true);
[foldiakObj, initWeights] = Foldiak(stimuliSets=foldiakParams, ...
                     weightsRange=0.1, ...
                     nComplex=4);
trainedWeights = foldiakObj.train(alpha=0.02, delta=0.2);

FoldiakGraphics.showConnections(8, 8, 4, 4, initWeights)
FoldiakGraphics.showConnections(8, 8, 4, 4, trainedWeights)
%% PLOT HISTOGRAMS OF WEIGHTS FOR EACH COMPLEX UNIT
weights = trainedWeights;
subplot(2,2,1)
bar(trainedWeights(:,1))
subplot(2,2,2)
bar(trainedWeights(:,2))
subplot(2,2,3)
bar(trainedWeights(:,3))
subplot(2,2,4)
bar(trainedWeights(:,4))
