time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
%% INITIALIZE DEFAULT DENDRITE MODEL - copied from hpoParallel.m
function modelNeuronObj = modelInit(seed, hps, hpPrev)
    % rng(seed, "twister"); % this needs to be reinitialized every time to keep results reproducible

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
seed = 0;
stepSize = 0.01; % fractional jump
gridSize = 100;
hps = {'scaleNMDA'};
hpInitVals = -0.1;
% hps = {'duPotent'; 'duDepress'; 'duDecay'; 'duBaseline'; 'scaleNMDA'; 'scaleNoNMDA'}; % hyperparameters to tune
% hpInitVals = [3 -0.3 -0.3 0 -0.1 0.003].';
nHPs = numel(hps);
%%
hpVals = hpInitVals * (1 + stepSize) .^ ((0:gridSize)-floor(gridSize/2));
lossVals = nan(nHPs, gridSize + 1);

for iHP = 1:nHPs
    tic
    fprintf("Grid search on %s (%d/%d)... \t", hps{iHP}, iHP, nHPs)
    parfor iGrid = 1:gridSize
        try
            hpPrev = hpInitVals;
            hpPrev(iHP) = hpVals(iHP, iGrid);
            modelNeuronObj = modelInit(seed, hps, hpPrev);
            [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
            lossVals(iHP, iGrid) = calcLoss(modelNeuronObj, resultsRFAfter);
        catch
            warning("Attempted to use hp %s = %1.3f", hps{iHP}, hpVals(iHP, iGrid))
        end
    end
    toc
end
%%
figure(Name="Grid search")
tl = tiledlayout(ceil(nHPs/2), 2);
for iHP = 1:nHPs
    nexttile
    % figure(Name=sprintf("%s", hps{iHP}));
    hold on
    plot(hpVals(iHP, :), lossVals(iHP, :))
    scatter(hpVals(iHP, :), lossVals(iHP, :), "filled", Color='red');
    % plot(hpValsPos(iHP, :), lossValsPos(iHP, :))
    % scatter(hpValsPos(iHP, :), lossValsPos(iHP, :), "filled", Color='blue');
    title(sprintf("%s", hps{iHP}));
    xlabel("HP value")
    ylabel("Loss")
    hold off
end
%% record values of hps and their loss
writematrix(hpVals, [time, '_hpVals','.txt'])
writematrix(lossVals, [time, '_lossVals','.txt'])
fclose('all');