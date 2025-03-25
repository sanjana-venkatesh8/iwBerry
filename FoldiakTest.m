% test script for Foldiak model
foldiakParams = StimuliParams(nX=8, nY=8, nOrient=4, barLength=8, scanLength=8, ...
    nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=400, ...
    nDendRecord=100, recordingTimeInterval=1, isFoldiak=true);
% stimuliParamsObj = StimuliParams(8, 8, 4, 8, 8, 0, 0, 0, 0, 1000, 400, 3, 1, true);
[foldiakObj, initWeights] = Foldiak(stimuliSets=foldiakParams, ...
                     weightsRange=0.1, ...
                     nComplex=100);
[trainedWeights, complexRF, cumulStimRF] = foldiakObj.train(alpha=0.02, delta=0.2);

% FoldiakGraphics.showConnections(foldiakParams.nX, foldiakParams.nY, foldiakParams.nOrient, foldiakObj.nComplex, initWeights)
% FoldiakGraphics.showConnections(foldiakParams.nX, foldiakParams.nY, foldiakParams.nOrient, foldiakObj.nComplex, trainedWeights)
%% PLOT HISTOGRAMS OF WEIGHTS FOR EACH COMPLEX UNIT
weights = trainedWeights;
for i = 1:foldiakObj.nComplex
    subplot(2, ceil(foldiakObj.nComplex/2), i)
    bar(trainedWeights(:, i))
    title(sprintf("Complex unit #%d", i))
end

%% ORIENTATION TUNING CURVES - NEED TO FIX
% calculate orientation tuning indices
complexUnitIOrient = zeros(foldiakObj.nComplex, 1);
compositeOrient = zeros(foldiakParams.nOrient, 1);
allBranchOrient = zeros(foldiakObj.nComplex, foldiakParams.nOrient);
branchIOrient = zeros(foldiakObj.nComplex, 1);

cumulOrient = zeros(foldiakParams.nOrient, 1);
for iOrient = 1:foldiakParams.nOrient
    for iX = 1:foldiakParams.nX
        for iY = 1:foldiakParams.nY
            cumulOrient(iOrient) = cumulOrient(iOrient) + ...
                cumulStimRF(iX, iY, iOrient);
        end
    end
end

for iComplex = 1:foldiakObj.nComplex
    orientTuning = zeros(foldiakParams.nOrient, 1);
    spatialRF = zeros(foldiakParams.nX, foldiakParams.nY);

    for iOrient = 1:foldiakParams.nOrient
        orientTuning(iOrient) = sum(complexRF(iComplex,:,:,iOrient), 'all');
    end
    allBranchOrient(iComplex, :) = orientTuning.' ./ cumulOrient.';
    compositeOrient = compositeOrient + orientTuning;
    % for iX = 1:foldiakParams.nX
    %     for iY = 1:foldiakParams.nY
    %         spatialRF(iX, iY) = sum(args.branchSpikeHist(iComplex, iX, iY, :), 'all');
    %     end
    % end
end

for iComplex = 1:foldiakObj.nComplex
    orientTuning = allBranchOrient(iComplex, :);
    pref = find(orientTuning == max(orientTuning), 1, 'first');
    null = find(orientTuning == min(orientTuning), 1, 'first');

    orientIndex = (orientTuning(pref) - orientTuning(null)) / ...
                    sum(orientTuning);
    branchIOrient(iComplex) = orientIndex;
end


% plot results
nTimesteps = foldiakParams.nTrainInst * foldiakParams.scanLength;
branchNMDASpikeRate = sum(complexRF, 2:4) / nTimesteps;

figure(Name="Layout 1: Orientation Tuning Across Complex Units (Foldiak model)")            
subplot(2, 1, 1)
scatter(branchNMDASpikeRate, ...
    branchIOrient, 50); % nG.resultsRFBefore.branchSize1(1:(end-2)), 'filled');
c = colorbar;
c.Label.String = "RF Size (pixels)";
title("Before plasticity")
ylabel("Orientation Tuning Index")
xlabel("NMDA Spike Rate")

% subplot(2, 1, 2)
% scatter(nG.resultsAfter.branchNMDASpikeRate, ...
%     nG.resultsRFAfter.branchIOrient(1:(end-2)), 50, nG.resultsRFAfter.branchSize1(1:(end-2)), 'filled');
% c = colorbar;
% c.Label.String = "RF Size (pixels)";
% title("After plasticity")
% ylabel("Orientation Tuning Index")
% xlabel("NMDA Spike Rate")
