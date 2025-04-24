%% DEFINE LOSS FUNCTION
% Loss (J) = c * 2 * (1 - orientationTuning) + (1 - c) * 1/64 * (64 - RFSize2)
% c = [0, 1] - weighting factor [CURRENTLY NOT IN USE]

function [totalLoss, perBranchLoss, somaLoss] = calcLoss(mN, resultsRFAfter)
    arguments
        mN ModelNeuron
        resultsRFAfter RFResults
    end

    % TODO: figure out if using branchSize1 or branchSize2
    % c = 0.5;
    fieldSize = mN.stimParams.nX * mN.stimParams.nY;
    perBranchLossIOrient = 1 - resultsRFAfter.branchIOrient(1:mN.dendParams.nBranches);
    perBranchLossRFSize = (fieldSize - resultsRFAfter.branchSize1(1:mN.dendParams.nBranches));
    
    perBranchLoss = 2 .* perBranchLossIOrient + (1/fieldSize) .* perBranchLossRFSize;
    perBranchLoss(isnan(perBranchLoss)) = 2; % fill in NaN values (i.e. branches that had 0 spikes) with 2 (maximum reasonable loss value)    
    totalLoss = sum(perBranchLoss, 'all');

    iSoma = length(resultsRFAfter.branchIOrient) - 1;
    somaLoss = 2 * resultsRFAfter.branchIOrient(iSoma) + (1/fieldSize) * resultsRFAfter.branchSize1(iSoma);
end