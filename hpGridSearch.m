%%
function loss = calcLoss(mN, resultsRFAfter)
    arguments
        mN ModelNeuron
        resultsRFAfter RFResults
    end

    % c = 0.5;
    fieldSize = mN.stimParams.nX * mN.stimParams.nY;
    loss = sum(1 - rmmissing(resultsRFAfter.branchIOrient), 'all') + ...
        1/(fieldSize) * sum(fieldSize - rmmissing(resultsRFAfter.branchSize1(1:mN.dendParams.nBranches)), 'all');
    % fprintf("DEBUG, Loss = %f\nDEBUG, # of branches with NaN tuning: %d.\n", loss, sum(isnan(resultsRFAfter.branchIOrient), 'all'))
end
%% INITIALIZE DEFAULT DENDRITE MODEL - copied from hpoParallel.m
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
%%
seed = 0;
stepSize = 0.01; % fractional jump
gridSize = 100;
hps = {'duPotent'; 'duDepress'; 'duDecay'; 'duBaseline'; 'scaleNMDA'; 'scaleNoNMDA'}; % hyperparameters to tune
hpInitVals = [3 -0.3 -0.3 0 -0.1 0.003].';
nHPs = numel(hps);

hpVals = hpInitVals * (1 + stepSize) .^ ((0:gridSize)-gridSize);
lossVals = nan(nHPs, gridSize + 1);

for iHP = 1:nHPs
    tic
    fprintf("Grid search on %s (%d/%d)... \t", hps{iHP}, iHP, nHPs)
    parfor iGrid = 1:gridSize
        hpPrev = hpInitVals;
        hpPrev(iHP) = hpVals(iHP, iGrid);
        modelNeuronObj = modelInit(seed, hps, hpPrev);
        [~, ~, ~, ~, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=false);
        lossVals(iHP, iGrid) = calcLoss(modelNeuronObj, resultsRFAfter);
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
    title(sprintf("%s", hps{iHP}));
    xlabel("HP value")
    ylabel("Loss")
    hold off
end
%% record values of hps and their loss
time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
writematrix(hpVals, [time, '_hpVals','.txt'])
writematrix(lossVals, [time, '_lossVals','.txt'])
fclose('all');