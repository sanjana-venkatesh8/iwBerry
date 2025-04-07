%% DEFINE LOSS FUNCTION
% Loss (J) = c * 2 * (1 - orientationTuning) + (1 - c) * 1/64 * (64 - RFSize2)
% c = [0, 1] - weighting factor [CURRENTLY NOT IN USE]

function loss = calcLoss(mN, resultsRFAfter)
    arguments
        mN ModelNeuron
        resultsRFAfter RFResults
    end

    % c = 0.5;
    fieldSize = mN.stimParams.nX * mN.stimParams.nY;
    lossIOrient = 1 - resultsRFAfter.branchIOrient(1:mN.dendParams.nBranches);
    lossRFSize = 1/(fieldSize) .* (fieldSize - resultsRFAfter.branchSize1(1:mN.dendParams.nBranches));
    loss = lossIOrient + lossRFSize;
    loss(isnan(loss)) = 2; % fill in NaN values (i.e. branches that had 0 spikes) with 2 (maximum reasonable loss value)    
    loss = sum(loss, 'all');
    % fprintf("DEBUG, Loss = %f\nDEBUG, # of branches with NaN tuning: %d.\n", loss, sum(isnan(resultsRFAfter.branchIOrient), 'all'))
end