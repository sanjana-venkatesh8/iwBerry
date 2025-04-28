for backpropAP = [0]
% INITIALIZATION
rng(0, "twister");
tic % start stopwatch to time model execution

dendParams = dendParamConfig.dendParamsBerry;
dendParams.backpropAP = backpropAP;
% dendParams.duPotent = 3.414042;
% dendParams.duDepress = -0.282876;
% dendParams.scaleNMDA = -0.129623;
% dendParams.scaleNoNMDA = 0.001866;

% CHANGED number of test/train instances - increase test size to 5k so loss is
% smoother, decrease train size to 2.5k for time efficiency
stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    noiseType='exact', noiseAmt=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

% SIMULATION
[resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=true);

toc % print elapsed time

fprintf("Loss = %1.3f\n", calcLoss(modelNeuronObj, resultsRFAfter));

% GRAPHICS
iBranch = 10; % randi(stimuliParamsObj.nDendRecord);
iSyn = 10; % randi(dendParams4Obj.branchSize);
tau = 150;
nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

% nG.layout1()
% nG.layout2(iBranch)
nG.layout3()
% nG.layout4(iBranch, dendParams.branchGLeak, tau)
% nG.layout5(iBranch, iSyn, dendParams.branchGLeak);
% nG.layout7()
end
%% 3/31/2025 noisy stimuli
rng(0, "twister");
dendParams = dendParamConfig.dendParamsBerry;

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    noiseType='exact', noiseAmt=1, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=1);

modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

nG.plotStimulus("iStim",1, "isTrain",true, "isPlotWithOrientation",true)