% Title: superlinear_bias_script.m
% Description: Calculate the bias in the branches as a function of the bias
% in stimuli (toward orientation 1)

time = char(datetime('now', 'Format', "MM-dd-yyyy-HHmm"));
nTrialReps = 1;
step_size = 0.005;
stim_bias = step_size*(0:25);
branch_bias = zeros(length(stim_bias), nTrialReps);

for i = stim_bias
    for rep = 1:nTrialReps
        % INITIALIZATION
        tic % start stopwatch to time model execution

        % rng(0, "twister"); % FREEZE NOISE
        dendParams = dendParamConfig.dendParamsBerry;
        dendParams.backpropAP = 0; % TURN OFF BAP
        dendParams.duPotent = 2;

        stimuliParamsObj = StimuliParams(nX=8, nY=8, nOrient=4, barLength=4, scanLength=4, ...
            noiseType='exact', noiseAmt=0, propNoiseScan=0, ...
            nTestRepeats=0, nTestInst=1000, nTrainInst=10000, ...
            nDendRecord=dendParams.nBranches, recordingTimeInterval=1, isFoldiak=false, isWraparound=1, ...
            orientPmf=[0.25+3*i, 0.25-i, 0.25-i, 0.25-i]);

        modelNeuronObj = ModelNeuron(dendParams=dendParams, stimParams=stimuliParamsObj, plasticityFlag=4);

        % SIMULATION
        [resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuronObj] = modelNeuronObj.jens2Plasticity(showOutput=true);

        branchPrefCounts = histcounts(resultsRFAfter.branchPref(1:end-2));
        branch_bias(int32(1/step_size*i)+1, rep) = branchPrefCounts(1)/dendParams.nBranches - 0.25;
    
        toc % print elapsed time
    end
end

writematrix(stim_bias, [time, '_stimBias','.txt'])
writematrix(branch_bias, [time, '_branchBias','.txt'])