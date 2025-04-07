% INITIALIZATION
rng(0, "twister");
tic % start stopwatch to time model execution

dendParams = dendParamConfig.dendParamsGridSearch;

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

% SIMULATION
[resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);

toc % print elapsed time

fprintf("Loss = %1.3f\n", calcLoss(modelNeuronObj, resultsRFAfter));

%% GRAPHICS
iBranch = 10; % randi(stimuliParamsObj.nDendRecord);
iSyn = 10; % randi(dendParams4Obj.branchSize);
tau = 150;
nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

nG.layout1()
nG.layout2(iBranch)
nG.layout3()
nG.layout4(iBranch, dendParams.branchGLeak, tau)
% nG.layout5(iBranch, iSyn, dendParams3Obj.branchGLeak);
nG.layout7()

%% 3/31/2025 noisy stimuli
rng(0, "twister");
dendParams = dendParamConfig.dendParamsBerry;

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    nL4ExactNoise=0, nL4PoissNoise=100, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=1);

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

nG.plotStimulus("iStim",1, "isTrain",true, "isPlotWithOrientation",true)
%% calcLoss()
% fieldSize = stimuliParamsObj.nX * stimuliParamsObj.nY;
% loss = 2 * sum(1 - resultsRFAfter.branchIOrient, 'all');% + ...
%     % 1/(fieldSize) * sum(fieldSize - resultsRFAfter.branchSize1(1:dendParams4Obj.nBranches), 'all');
% loss = loss / dendParams.nBranches;
% fprintf("Loss = %f\n", loss);

%% LOSS FUNCTION (2.17.2025)
% bad performance (8/64 squares covered) - J = 
% --> ~0.1 (orientationTuning)
% --> ~ 0.875 (RFSize2)
% good performance (24/64 squares covered) - J = 
% --> 0 (orientationTuning)
% --> ~ 0.6 (RFSize2)

% Loss (J) = c * 2(?) * (1 - orientationTuning) + (1 - c) * 1/64 * (64 - RFSize2)
% c = [0, 1] - weighting factor
% function loss = calcLoss()
%     c = 0.5;
%     fieldSize = stimuliParamsObj.nX * stimuliParamsObj.nY;
%     loss = sum(1 - resultsRFAfter.branchIOrient, 'all') + ...
%         1/(fieldSize) * sum(fieldSize - resultsRFAfter.branchSize1(1:dendParams4Obj.nBranches), 'all');
% end