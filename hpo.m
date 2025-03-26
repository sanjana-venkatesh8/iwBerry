% Author: Sanjana Venkatesh <sv8281@princeton.edu>
% Date: March 2025
% Description: 
%   Gradient descent optimization of hyperparameters to the Jens
%   dendrite model. Tune branch parameters before somatic parameters.
%   Hyperparameters:
%       - branch parameters:  duPotent, duDepress, duDecay, duBaseline, scaleNMDA, scaleNoNMDA
%       - somatic parameters: backpropAP

nGDSteps = 40; % number of gradient descent steps
a = 0.01; % gradient descent learning rate
b = 0.01; % fractional size of hyperparameter perturbation when calculating gradient

hps = {'duPotent'; 'duDepress'; 'duDecay'; 'duBaseline'; 'scaleNMDA'; 'scaleNoNMDA'}; % hyperparameters to tune
%% INITIALIZE DEFAULT DENDRITE MODEL
function modelNeuronObj = modelInit(seed)
    rng(seed, "twister"); % this needs to be reinitialized every time to keep results reproducible

    dendParamsDefault = DendriteParams(nSynapses=2000, branchSize=20, ...
        synStrength=10, synAttenuation=10, ...
        branchThresh=15, branchStrength=7, somaThresh=20, backpropAP=7, ...
        synNoise=0.2, branchNoise=0.2, somaNoise=0, weightsRange=1, ...
        duPotent=3, duDepress=-0.3, duDecay=-0.3, duBaseline=0, ...
        uRecycle=-3, uMax=10, plastTime=4, ...
        scaleNMDA=-0.1, scaleNoNMDA=0.003, branchGLeak=1, ...
        initSynGInhib=4, somaPlus=0.1, somaMinus=-0.005, somaGLeak=0.1, initSomaGInhib=0.5);
    
    stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
        nL4ExactNoise=0, nL4PoissNoise=0, propNoiseScan=0, ...
        nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
        nDendRecord=dendParamsDefault.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);
    
    modelNeuronObj = ModelNeuron(dendParams=dendParamsDefault, stimParams=stimuliParamsObj, plasticityFlag=4);
end
%% DEFINE LOSS CALCULATION FUNCTION
% Loss (J) = c * 2 * (1 - orientationTuning) + (1 - c) * 1/64 * (64 - RFSize2)
% c = [0, 1] - weighting factor [CURRENTLY NOT IN USE]

function loss = calcLoss(mN, resultsRFAfter)
    arguments
        mN ModelNeuron
        resultsRFAfter RFResults
    end

    % c = 0.5;
    fieldSize = mN.stimParams.nX * mN.stimParams.nY;
    loss = sum(1 - rmmissing(resultsRFAfter.branchIOrient), 'all') + ...
        1/(fieldSize) * sum(fieldSize - rmmissing(resultsRFAfter.branchSize1(1:mN.dendParams.nBranches)), 'all');
    fprintf("DEBUG, Loss = %f\nDEBUG, # of branches with NaN tuning: %d.\n", loss, sum(isnan(resultsRFAfter.branchIOrient), 'all'))
end
%% NOISE TESTING
JRecord = NaN(nGDSteps ,1); % value of loss (J) at each GD step

% create list of random seeds - DEBUG: freeze noise completely during GD
rng('shuffle', "twister")
seedlist = randi(100000,nGDSteps,1);

for iStep = 1:nGDSteps
    tic
    
    % initialize model
    modelNeuronObj = modelInit(seedlist(iStep));

    % run model with plasticity and store loss
    [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
    J0 = calcLoss(modelNeuronObj, resultsRFAfter);
    JRecord(iStep) = J0;

    toc
end

time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
writematrix(JRecord, [time, 'JRecord_','.txt']);
%% GRADIENT DESCENT
% JRecord = NaN(nGDSteps ,1); % value of loss (J) at each GD step
% JGradRecord = NaN(nGDSteps, size(hps, 1)); % value of loss gradient at each GD step
% hpRecord = NaN(nGDSteps + 1, size(hps, 1)); % value of each hyperparameter at each GD step
% 
% % create file to save parameter values
% time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
% % hpFileID = fopen(['hps_', time,'.txt'], 'w');
% 
% % create list of random seeds - DEBUG: freeze noise completely during GD
% rng('shuffle', "twister")
% seedlist = zeros(nGDSteps, 1); % randi(100000,nGDSteps,1);
% 
% for iStep = 1:nGDSteps
%     tic
% 
%     % initialize model
%     modelNeuronObj = modelInit(seedlist(iStep));
% 
%     % store initial HP values
%     if iStep == 1
%         for iHP = 1:size(hps, 1)
%             hpRecord(1, iHP) = modelNeuronObj.dendParams.(hps{iHP});
%         end
%     end
% 
%     % sub in previous step's HP values to model
%     for iHP = 1:size(hps, 1)
%         modelNeuronObj.dendParams.(hps{iHP}) = hpRecord(iStep, iHP);
%     end
%     % run model with plasticity and store loss
%     [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
%     J0 = calcLoss(modelNeuronObj, resultsRFAfter);
%     JRecord(iStep) = J0;
% 
%     % numerically calculate gradient by perturbing each of the hyperparameters
%     for iHP = 1:size(hps, 1)
% 
%         fprintf("\tGD Step %d, perturbing modelNeuronObj.dendParams.%s\n", iStep, hps{iHP})
% 
%         % reinitialize model, increase the hyperparameter by a percentage b
%         modelNeuronObj = modelInit(seedlist(iStep));
%         % sub in previous step's HP values to model
%         for jHP = 1:size(hps, 1)
%             modelNeuronObj.dendParams.(hps{jHP}) = hpRecord(iStep, jHP);
%         end
%         % perturb the HP
%         modelNeuronObj.dendParams.(hps{iHP}) = (1 + b) * hpRecord(iStep, iHP);
% 
%         % run model and calculate loss
%         [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
%         Ji = calcLoss(modelNeuronObj, resultsRFAfter);
% 
%         % store value of partial derivative with respect to this
%         % hyperparameter
%         JGradRecord(iStep, iHP) = Ji - J0;
%     end
% 
%     % GD update: hp_(i+1) <= hp_i - a * grad(hp_i). Store HP values at this
%     % step
%     for iHP = 1:size(hps, 1)
%         hpRecord(iStep + 1, iHP) = ...
%             modelNeuronObj.dendParams.(hps{iHP}) - a * JGradRecord(iStep, iHP);
%     end
% 
%     fprintf("GD Step %d: Loss = %1.3f, ", iStep, J0)
%     toc
% end
% 
% %% record values of loss and its gradient per GD step
% writematrix(JRecord, [time, 'JRecord_','.txt'])
% writematrix(JGradRecord, [time, 'JGradRecord_','.txt'])
% writematrix(hpRecord, [time, 'hps_','.txt'])
% fclose('all');