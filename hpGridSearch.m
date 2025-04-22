time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
%%
isNoiseFrozen = false;
seed = 0;
stepSize = 0.01; % fractional jump
gridSize = 200;
hps = {'scaleNMDA'};
hpInitVals = -0.1;
    % hps = {'duPotent'; 'duDepress'; 'duDecay'; 'duBaseline'; 'scaleNMDA'; 'scaleNoNMDA'}; % hyperparameters to tune
    % hpInitVals = [3 -0.3 -0.3 0 -0.1 0.003].';
nHPs = numel(hps);
nTrialReps = 1;
%% INITIALIZE DEFAULT DENDRITE MODEL
function modelNeuronObj = modelInit(isNoiseFrozen, seed, hps, hpPrev)
    if isNoiseFrozen
        rng(seed, "twister"); % this needs to be reinitialized every time to keep results reproducible
    end

    dendParams = dendParamConfig.dendParamsBerry;

    % CHANGED number of test/train instances - increase test size to 5k so loss is
    % smoother, decrease train size to 2.5k for time efficiency
    stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
        noiseType='exact', noiseAmt=0, propNoiseScan=0, ...
        nTestRepeats=0, nTestInst=5000, nTrainInst=2500, ...
        nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=0);
    
    modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

    if exist('hps', 'var') && exist('hpPrev', 'var')
        for iHP = 1:size(hps, 1)
            modelNeuronObj.dendParams.(hps{iHP}) = hpPrev(iHP);
        end
    end
end
%%
function [lossAvg, lossStdDev] = calcLossStats(nTrialReps, isNoiseFrozen, seed, hps, hpPrev)
    repLoss = nan(nTrialReps, 1);
    parfor iRep = 1:nTrialReps % run each trial several times and average the loss
        modelNeuronObj = modelInit(isNoiseFrozen, seed, hps, hpPrev);
        [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
        repLoss(iRep) = calcLoss(modelNeuronObj, resultsRFAfter);
    end

    lossAvg = mean(repLoss);
    lossStdDev = std(repLoss);
end
%%
hpVals = hpInitVals * (1 + stepSize) .^ ((0:gridSize)-floor(gridSize/2));
lossVals = nan(nHPs, gridSize + 1);
lossStdDevVals = nan(nHPs, gridSize + 1);

for iHP = 1:nHPs
    tic
    fprintf("Grid search on %s (%d/%d)... \t", hps{iHP}, iHP, nHPs)
    parfor iGrid = 1:gridSize
        try
            hpPrev = hpInitVals;
            hpPrev(iHP) = hpVals(iHP, iGrid);
            [lossVals(iHP, iGrid), lossStdDevVals(iHP, iGrid)] = calcLossStats(nTrialReps, isNoiseFrozen, seed, hps, hpPrev);
        catch
            warning("Attempted to use hp %s = %1.3f", hps{iHP}, hpVals(iHP, iGrid))
        end
    end
    toc
end
%%
figure(Name="Grid search")
% tl = tiledlayout(ceil(nHPs/2), 2);
for iHP = 1:nHPs
    nexttile
    % figure(Name=sprintf("%s", hps{iHP}));
    hold on
    plot(hpVals(iHP, :), lossVals(iHP, :))
    scatter(hpVals(iHP, :), lossVals(iHP, :), "filled", Color='red');
    if nTrialReps > 1
        errorbar(hpVals(iHP, :), lossVals(iHP, :), lossStdDevVals(iHP, :))
    end
    title(sprintf("%s", hps{iHP}));
    xlabel("HP value")
    ylabel("Loss")
    hold off
end
%% record values of hps and their loss
writematrix(hpVals, [time, '_hpVals','.txt'])
writematrix(lossVals, [time, '_lossVals','.txt'])
if nTrialReps > 1
    writematrix(lossStdDevVals, [time, '_lossStdDevVals', '.txt'])
end
fclose('all');