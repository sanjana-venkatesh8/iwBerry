rng(0, "twister");
tic % start stopwatch to time model execution

dendParamsDefault = DendriteParams(nSynapses=2000, branchSize=20, ...
    synStrength=10, synAttenuation=10, ...
    branchThresh=15, branchStrength=7, somaThresh=20, backpropAP=7, ...
    synNoise=0.2, branchNoise=0.2, somaNoise=0, weightsRange=1, ...
    duPotent=3, duDepress=-0.3, duDecay=-0.3, duBaseline=0, ...
    uRecycle=-3, uMax=10, plastTime=4, ...
    scaleNMDA=-0.1, scaleNoNMDA=0.003, branchGLeak=1, ...
    initSynGInhib=4, somaPlus=0.1, somaMinus=-0.005, somaGLeak=0.1, initSomaGInhib=0.5);

% dendParamsDefault = DendriteParams(nSynapses=2000, branchSize=20, ...
%     synStrength=10, synAttenuation=10, ...
%     branchThresh=15, branchStrength=7, somaThresh=20, backpropAP=7, ...
%     synNoise=0.2, branchNoise=0.2, somaNoise=0, weightsRange=1, ...
%     duPotent=3.011, duDepress=-0.3049, duDecay=-0.3084, duBaseline=0, ...
%     uRecycle=-3, uMax=10, plastTime=4, ...
%     scaleNMDA=-0.1072, scaleNoNMDA=0.0091, branchGLeak=1, ...
%     initSynGInhib=4, somaPlus=0.1, somaMinus=-0.005, somaGLeak=0.1, initSomaGInhib=0.5);

stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
    nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
    nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
    nDendRecord=dendParamsDefault.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);

modelNeuronObj = ModelNeuron(dendParams=dendParamsDefault, stimParams=stimuliParamsObj, plasticityFlag=4);

[resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);

toc % print elapsed time

%% GRAPHICS
iBranch = 10; % randi(stimuliParamsObj.nDendRecord);
iSyn = 10; % randi(dendParams4Obj.branchSize);
tau = 150;
nG = NeuronGraphics(modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter);

nG.layout1()
nG.layout2(iBranch)
nG.layout3()
nG.layout4(iBranch, dendParamsDefault.branchGLeak, tau)
% nG.layout5(iBranch, iSyn, dendParams3Obj.branchGLeak);
nG.layout7()

% calcLoss()
fieldSize = stimuliParamsObj.nX * stimuliParamsObj.nY;
loss = 2 * sum(1 - resultsRFAfter.branchIOrient, 'all');% + ...
    % 1/(fieldSize) * sum(fieldSize - resultsRFAfter.branchSize1(1:dendParams4Obj.nBranches), 'all');
loss = loss / dendParamsDefault.nBranches;
fprintf("Loss = %f\n", loss);

%% LOSS FUNCTION (2.17.2025)
% bad performance (8/64 squares covered) - J = 
% --> ~0.1 (orientationTuning)
% --> ~ 0.875 (RFSize2)
% good performance (24/64 squares covered) - J = 
% --> 0 (orientationTuning)
% --> ~ 0.6 (RFSize2)

% Loss (J) = c * 2(?) * (1 - orientationTuning) + (1 - c) * 1/64 * (64 - RFSize2)
% c = [0, 1] - weighting factor
function loss = calcLoss()
    c = 0.5;
    fieldSize = stimuliParamsObj.nX * stimuliParamsObj.nY;
    loss = sum(1 - resultsRFAfter.branchIOrient, 'all') + ...
        1/(fieldSize) * sum(fieldSize - resultsRFAfter.branchSize1(1:dendParams4Obj.nBranches), 'all');
end