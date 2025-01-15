% Make a scan and display the L4 activity/spatial location of the bar at 
% t = 0.
stimuliParamsObj = StimuliParams(8, 8, 4, 4, 4, 0, 0, 0, 0, 5, 10, 3, 1)
%% TODO: TURN THIS INTO A GRAPHICS FUNCTION
[L4Activity, barXLoc, barYLoc] = makeScan(stimParams=stimuliParamsObj, showOutput=true);

for iOrient = 1:stimuliParamsObj.nOrient
    L4ChooseOrient = L4Activity(iOrient:4:256,1);
    L4Heat = reshape(L4ChooseOrient, 8, 8).';
    subplot(2, 2, iOrient)
    title(Name=sprintf("L4 neurons with orientation %d at t = 0", iOrient))
    heatmap(L4Heat())
    % colormap('hot')
end
 
for iScan = 1
    figure
    h1 = axes;
        scatter(barXLoc(:,iScan), barYLoc(:,iScan), "filled")
    set(h1,'YDir','reverse')
    xlim([1 8])
    ylim([1 8])
end
%%
tic % start stopwatch to time model execution

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=2500, ...
    nDendRecord=3, recordingTimeInterval=1, isFoldiak=false);
% stimuliParamsObj = StimuliParams(8, 8, 4, 4, 4, 0, 0, 0, 0, 1000, 10000, 3, 1);
% dendParamsObj = DendriteParams(2000, 10, 20, 10, 10, 20, 7, 20, 0, 0.2, 0.2, 0, 1, 2, -0.2, -0.1, 0, -3, 10, 4, -0.1, 0.003, 0, 1, 1, -0.05, 0, 1);
% without noise:
% dendParamsObj = DendriteParams(2000, 10, 20, 10, 10, 20, 7, 20, 0, 0, 0, 0, 1, 2, -0.2, -0.1, 0, -3, 10, 4, -0.1, 0.003, 0, 1, 1, -0.05, 0, 1);
dendParamsObj = DendriteParams(2000, 10, 20, 10, 10, 15, 7, 20, 0, 0.2, 0.2, 0, 1, 2, -0.3, -0.3, 0, -3, 10, 4, -0.1, 0.003, 1, 4, 0.1, -0.005, 0.1, 0.5);
% dendParamsObj = DendriteParams(2000, 10, 20, 10, 10, 15, 7, 20, 0, 0.2, 0.2, 0, 1, 2, -0.3, -0.3, 0, -3, 10, 4, -0.1, 0.003, 1, 4, 0.1, -0.005, 0.1, 0.5);
modelNeuronObj = ModelNeuron(dendParams=dendParamsObj, stimParams=stimuliParamsObj, plasticityFlag=4);
[resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity();

toc % print elapsed time
%% GRAPHICS
iBranch = randi(stimuliParamsObj.nDendRecord);
iSyn = randi(dendParamsObj.branchSize);
tau = 150;
nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);
% nG.layout1()
% nG.layout2(iBranch)
nG.layout4(iBranch, dendParamsObj.branchGLeak, tau)
nG.layout5(iBranch, iSyn, dendParamsObj.branchGLeak);
%%
[allBranchActivity, allBranchGExc, allBranchGain, ...
    allBranchDidSpike, totalNMDASpikes, synActivity, modelNeuron] = ...
    modelNeuronObj.calcDendActivity(plasticityFlag=1, doBackprop=true, stimBarL4Activity=stimuliParamsObj.trainSet.L4Activity(:, 1));
%% MODEL TESTING
% create ModelNeuron object
% run the runModel method
%% SCRATCH