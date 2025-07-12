% INITIALIZATION
tic % start stopwatch to time model execution

rng(37, "twister"); % FREEZE NOISE
dendParams = dendParamConfig.dendParamsBerry;
dendParams.backpropAP = 0; % TURN OFF BAP
dendParams.duPotent = 2;
dendParams.branchGLeak = 2;
dendParams.duDepress = -0.4;
% dendParams.duDecay = -0.1;
% dendParams.scaleNMDA = -0.15;
% dendParams.somaThresh = 15;

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    noiseType='exact', noiseAmt=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=1, ...
    orientPmf='rand');

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

% SIMULATION
[resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=true);

toc % print elapsed time

% fprintf("Loss = %1.3f\n", calcLoss(modelNeuronObj, resultsRFAfter));
[totalLoss, branchLoss, somaLoss] = calcLoss(modelNeuronObj, resultsRFAfter);

% GRAPHICS
iBranch = 70; % randi(stimuliParamsObj.nDendRecord);
iSyn = 10; % randi(dendParams4Obj.branchSize);
nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);
nG.branchPrefHist();
nG.branchLockIn(iBranch=iBranch)

% nG.layout1()
% nG.layout2(iBranch)
nG.layout3(normalize=false)
% nG.layout4(iBranch, dendParams.branchGLeak)
% % nG.layout5(iBranch, iSyn, dendParams.branchGLeak);
% nG.layout7()
%% 3/31/2025 noisy stimuli
rng(0, "twister");
dendParams = dendParamConfig.dendParamsBerry;

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    noiseType='exact', noiseAmt=1, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=false);

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

nG.plotStimulus("iStim",14, "isTrain",true, "isPlotWithOrientation",true)