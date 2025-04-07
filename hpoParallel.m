% Author: Sanjana Venkatesh <sv8281@princeton.edu>
% Date: March 2025
% Description: 
%   PARALLELIZED gradient descent optimization of hyperparameters to the Jens
%   dendrite model. Tune branch parameters before somatic parameters.
%   Hyperparameters:
%       - branch parameters:  duPotent, duDepress, duDecay, duBaseline, scaleNMDA, scaleNoNMDA
%       - somatic parameters: backpropAP

nGDSteps = 40; % number of gradient descent steps
a = 10^-9; % gradient descent learning rate
b = 10^-5;%5 * 10^-6; % fractional size of hyperparameter perturbation when calculating gradient

% hps = {'duPotent'; 'duDepress'; 'duDecay'; 'duBaseline'; 'scaleNMDA'; 'scaleNoNMDA'}; % hyperparameters to tune
hps = {'duDepress'};%; 'duDepress'; 'duDecay'}; % DEBUG
%% INITIALIZE DEFAULT DENDRITE MODEL
function modelNeuronObj = modelInit(seed, hps, hpPrev)
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
        nTestRepeats=0, nTestInst=2500, nTrainInst=2500, ...
        nDendRecord=dendParamsDefault.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);
    
    modelNeuronObj = ModelNeuron(dendParams=dendParamsDefault, stimParams=stimuliParamsObj, plasticityFlag=4);

    if exist('hps', 'var') && exist('hpPrev', 'var')
        for iHP = 1:size(hps, 1)
            modelNeuronObj.dendParams.(hps{iHP}) = hpPrev(iHP);
        end
    end
end
% %% DEFINE LOSS FUNCTION
% % Loss (J) = 2 * (1 - orientationTuning) + 1/64 * (64 - RFSize2)
% % c = [0, 1] - weighting factor [CURRENTLY NOT IN USE]
% 
% function loss = calcLoss(mN, resultsRFAfter)
%     arguments
%         mN ModelNeuron
%         resultsRFAfter RFResults
%     end
% 
%     
%     fieldSize = mN.stimParams.nX * mN.stimParams.nY;
%     loss = sum(1 - rmmissing(resultsRFAfter.branchIOrient), 'all') + ...
%         1/(fieldSize) * sum(fieldSize - rmmissing(resultsRFAfter.branchSize1(1:mN.dendParams.nBranches)), 'all');
%     % fprintf("DEBUG, Loss = %f\nDEBUG, # of branches with NaN tuning: %d.\n", loss, sum(isnan(resultsRFAfter.branchIOrient), 'all'))
% end
%% GRADIENT DESCENT
JRecord = NaN(nGDSteps ,1); % value of loss (J) at each GD step
JGradRecord = NaN(nGDSteps, size(hps, 1)); % value of loss gradient at each GD step
hpRecord = NaN(nGDSteps + 1, size(hps, 1)); % value of each hyperparameter at each GD step

% create file to save parameter values
time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
% hpFileID = fopen(['hps_', time,'.txt'], 'w');

% create list of random seeds - DEBUG: freeze noise completely during GD
rng('shuffle', "twister")
seedlist = zeros(nGDSteps, 1); % randi(100000,nGDSteps,1);

for iStep = 1:nGDSteps
    tic

    seed = seedlist(iStep);

    % initialize model and store initial HP values if needed
    if iStep == 1
        modelNeuronObj = modelInit(seed);
        for iHP = 1:size(hps, 1)
            hpRecord(iStep, iHP) = modelNeuronObj.dendParams.(hps{iHP});
        end
        hpPrev = hpRecord(iStep, :);
    else
        hpPrev = hpRecord(iStep, :);
        modelNeuronObj = modelInit(seed, hps, hpPrev);
    end

    % run model with plasticity and store loss
    fprintf("\tGD Step %d, Calculating loss J0", iStep)
    [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
    J0 = calcLoss(modelNeuronObj, resultsRFAfter);
    JRecord(iStep) = J0;
    fprintf(" = %1.3f\n", J0)

    % numerically calculate gradient by perturbing each of the hyperparameters
    parfor iHP = 1:size(hps, 1)

        fprintf("\tGD Step %d, perturbing modelNeuronObj.dendParams.%s\n", iStep, hps{iHP})

        % reinitialize model, increase the hyperparameter by a percentage b
        modelNeuronObj = modelInit(seed, hps, hpPrev);

        % perturb the HP
        try
            modelNeuronObj.dendParams.(hps{iHP}) = (1 + b) * hpRecord(iStep, iHP);
            delta = abs(b * hpRecord(iStep, iHP)); % always a positive perturbation, but HP value might be negative
        catch
            warning("Tried to make %s = %1.3f [illegal]. Leaving value unchanged.", hps{iHP}, (1 + b) * hpRecord(iStep, iHP));
            delta = 1; % DEBUG - delta is not set.
        end

        % run model and calculate loss
        [resultsBefore, ~, resultsAfter, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
        Ji = calcLoss(modelNeuronObj, resultsRFAfter);
        spikeDiff = abs(sum(resultsAfter.NMDASpikesPerTimestep - resultsBefore.NMDASpikesPerTimestep, "all"));
        totalBits = size(resultsAfter.NMDASpikesPerTimestep, 1) * modelNeuronObj.dendParams.nBranches;
        fprintf("\tGD Step %d, # of bits flipped: %d/%d = %1.3f %% \n", iStep, spikeDiff, totalBits, 100*spikeDiff/totalBits)

        % store value of partial derivative with respect to this
        % hyperparameter
        % goes to inf
        JGradRecord(iStep, iHP) = (Ji - J0) / delta;
        % DEBUG
        fprintf("Ji %s = %1.3f. Jgrad = %1.3f\n", hps{iHP}, Ji, JGradRecord(iStep, iHP));
    end

    % GD update: hp_(i+1) <= hp_i - a * grad(hp_i). Store HP values at this
    % step. If they are illegal, leave the HP value unchanged.
    hpNext = hpPrev;
    parfor iHP = 1:size(hps, 1)
        modelNeuronObj = modelInit(seed, hps, hpPrev);
        try
            updatedVal = modelNeuronObj.dendParams.(hps{iHP}) - a * JGradRecord(iStep, iHP);
            modelNeuronObj.dendParams.(hps{iHP}) = updatedVal;
            hpNext(iHP) = updatedVal;
            % hpRecord(iStep + 1, iHP) = updatedVal;
        catch
            warning("Tried to make %s = %1.3f [illegal]. Leaving value unchanged", hps{iHP}, updatedVal);
            % hpRecord(iStep + 1, iHP) = modelNeuronObj.dendParams.(hps{iHP});
        end
    end
    hpRecord(iStep + 1, :) = hpNext;

    fprintf("GD Step %d: ", iStep)
    toc
end

%% record values of loss and its gradient per GD step
writematrix(JRecord, [time, '_JRecord_','.txt'])
writematrix(JGradRecord, [time, '_JGradRecord_','.txt'])
writematrix(hpRecord, [time, '_hps_','.txt'])
fclose('all');