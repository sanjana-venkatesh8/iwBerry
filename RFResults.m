classdef RFResults
    %RFResults Stores results after running dendriticRFAnalyze
    
    properties
        allBranchSpatial
        allBranchOrient
        branchIOrient % orient tuning index for each branch
        branchPref
        branchMaxResp
        branchVector
        branchSize1
        branchSize2
    end
    
    methods
        function RFResultsObj = RFResults(allBranchSpatial, ...
                allBranchOrient, branchIOrient, branchPref, ...
                branchMaxResp, branchVector, branchSize1, branchSize2)
            %RFResults Construct an instance of this class
            RFResultsObj.allBranchSpatial = allBranchSpatial;
            RFResultsObj.allBranchOrient = allBranchOrient;
            RFResultsObj.branchIOrient = branchIOrient;
            RFResultsObj.branchPref = branchPref;
            RFResultsObj.branchMaxResp = branchMaxResp;
            RFResultsObj.branchVector = branchVector;
            RFResultsObj.branchSize1 = branchSize1;
            RFResultsObj.branchSize2 = branchSize2;
        end
    end
end

