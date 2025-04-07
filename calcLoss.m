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
    loss = sum(1 - rmmissing(resultsRFAfter.branchIOrient), 'all') + ...
        1/(fieldSize) * sum(fieldSize - rmmissing(resultsRFAfter.branchSize1(1:mN.dendParams.nBranches)), 'all');
    % fprintf("DEBUG, Loss = %f\nDEBUG, # of branches with NaN tuning: %d.\n", loss, sum(isnan(resultsRFAfter.branchIOrient), 'all'))
end