classdef StimuliParams
    %STIMULIPARAMS Stores parameters to generate a set of stimuli
    % Stimulus parameters:
    %     nX = (int) # of x locations
    %     nY = (int) # of y locations
    %     nOrient = (int) # of possible bar orientations
    %     barLength = (int) bar length
    %     scanLength = (int) scan length
    %     nL4ExactNoise
    %     nL4PoissNoise
    %     propNoiseScan = (float, [0,1]) proportion of noise-only scans
    %     nTestRepeats
    %     nTestInst = (int) # of bars in test set
    %     nTrainInst = (int) # of bars in train set
    %     nDendRecord
    %     recordingTimeInterval
    %     trainSet
    %     testSet
    
    properties
        % TODO: datatype validation
        nX %int32
        nY %int32
        nOrient {mustBeMember(nOrient, [1 2 3 4])} %int32 % model only supports 0˚, 45˚, 90˚, 135˚
        barLength %int32
        scanLength %int32;
        noiseType;
        noiseAmt;
        propNoiseScan %double {mustBeInRange(propNoiseScan, 0, 1)}
        nTestRepeats %int32
        nTestInst %int32
        nTrainInst %int32
        nDendRecord %int32
        recordingTimeInterval;
        nL4Inputs;
        isFoldiak;
        isWraparound;
        % non-input properties
        trainSet;
        testSet;
    end
    
    methods (Access = public)
        function stimParams = StimuliParams(args)
            % nX, nY, ...
            %     nOrient, barLength, scanLength, nL4ExactNoise, ...
            %     nL4PoissNoise, propNoiseScan, nTestRepeats, nTestInst, ...
            %     nTrainInst, nDendRecord, recordingTimeInterval, isFoldiak, isWraparound)
            
            %STIMULIPARAMS Construct an instance of this class

            arguments
                args.nX; args.nY; args.nOrient; args.barLength; args.scanLength;
                args.noiseType; args.noiseAmt; args.propNoiseScan;
                args.nTestRepeats; args.nTestInst; args.nTrainInst;
                args.nDendRecord; args.recordingTimeInterval;
                args.isFoldiak = false; args.isWraparound = false;
            end

            try
                stimParams.nX = args.nX;
                stimParams.nY = args.nY;
                stimParams.nOrient = args.nOrient;
                stimParams.barLength = args.barLength;
                stimParams.scanLength = args.scanLength;
    
                stimParams.noiseType = validatestring(args.noiseType, {'exact', 'poisson'});
                stimParams.noiseAmt = args.noiseAmt;
                stimParams.propNoiseScan = args.propNoiseScan;

                stimParams.nTestRepeats = args.nTestRepeats; % variable not used for anything
                stimParams.nTestInst = args.nTestInst;
                stimParams.nTrainInst = args.nTrainInst;

                stimParams.nDendRecord = args.nDendRecord;
                stimParams.recordingTimeInterval = args.recordingTimeInterval;

                stimParams.isFoldiak = args.isFoldiak;
                stimParams.isWraparound = args.isWraparound;
                % % check if Foldiak model is creating this object. Foldiak
                % % model follows different indexing for L4/simple cells.
                % if ~exist(,'var')
                %    stimParams.isFoldiak = false;
                % else
                %    stimParams.isFoldiak = args.isFoldiak;
                % end

            catch
                throw(MException('StimuliParams:VarNotFound', ...
                    'One or more variables not found'))
            end

            % TODO: include showOutput param
            stimParams.nL4Inputs = args.nX * args.nY * args.nOrient;
            [stimParams.trainSet, stimParams.testSet] = makeTrainTestSets(stimParams=stimParams, showOutput=false);
        end
    % TODO: MAKE PRIVATE methods
    % end
    % 
    % methods(Access=private)
        % TODO: insert functions here if needed

        % Generate training and testing sets of scanning bar stimuli
        % Returns train and test structs each containing L4 activity, X and
        % Y locations for each scan in the set
        % TODO: make private method
        function [trainSet, testSet] = makeTrainTestSets(args) % equivalent of Do_Bar_TestSet
            arguments
                args.stimParams StimuliParams
                args.showOutput logical = false
            end

            stimParams = args.stimParams;
            showOutput = args.showOutput;

            trainSet.L4Activity = zeros(stimParams.nL4Inputs, ...
                stimParams.scanLength, stimParams.nTrainInst);
            trainSet.barXLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTrainInst); 
            trainSet.barYLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTrainInst); 
            trainSet.orientation = zeros(stimParams.nTrainInst, 1);
            
            testSet.L4Activity = zeros(stimParams.nL4Inputs, ...
                stimParams.scanLength, stimParams.nTestInst);
            testSet.barXLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTestInst); 
            testSet.barYLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTestInst); 
            testSet.orientation = zeros(stimParams.nTestInst, 1);

            % generate training set
            for iTrainInst = 1:stimParams.nTrainInst
                [trainSet.L4Activity(:, :, iTrainInst), ...
                    trainSet.barXLoc(:, :, iTrainInst), ...
                    trainSet.barYLoc(:, :, iTrainInst), ...
                    trainSet.orientation(iTrainInst)] ...
                                        = makeScan(stimParams=stimParams, ...
                                                   showOutput=showOutput);
            end

            % generate testing set
            for iTestInst = 1:stimParams.nTestInst
                [testSet.L4Activity(:, :, iTestInst), ...
                    testSet.barXLoc(:, :, iTestInst), ...
                    testSet.barYLoc(:, :, iTestInst), ...
                    testSet.orientation(iTestInst)] ...
                                        = makeScan(stimParams=stimParams, ...
                                                   showOutput=showOutput);
            end

        end

        % generate one scan and its corresponding L4 activity and X/Y
        % locations at each timestep
        function [L4Activity, barXLoc, barYLoc, orientation] = makeScan(args) % equivalent of Do_Bar_Scan
            arguments
                args.stimParams StimuliParams
                args.showOutput logical = false
            end

            stimParams = args.stimParams;
            showOutput = args.showOutput;

            L4Activity = zeros(stimParams.nL4Inputs, ...
                stimParams.scanLength); % stores which inputs are stimulated at each step of scan
            barXLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength); % stores x location of each point in bar at each step of scan
            barYLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength); % stores y location of each point in bar at each step of scan 

            isNoiseImage = (rand() < stimParams.propNoiseScan); % DEBUG 3/31: add this in
            if not(isNoiseImage)
                direction = 2 * randi(2) - 3; % set random scan direction (-1 or 1)
                orientation = randi(stimParams.nOrient); % set random orientation of bar ([1, nOrient])

                if showOutput
                    fprintf("Orientation = %d, direction = %d\n", orientation, direction)
                end

                % set random initial bar location
                if stimParams.isWraparound
                    xInit = randi(stimParams.nX);
                    yInit = randi(stimParams.nY);
                else
                    % start at a left/up offset to allow for truncations 
                    % on left/upper side as well as the right/lower sides 
                    % of the visual field.
                    
                    offset = 0;
                    % offset = stimParams.barLength - 2; % TODO: find out
                    % whether or not to include offset
                    xInit = randi(stimParams.nX + offset) - offset;
                    yInit = randi(stimParams.nY + offset) - offset;
                end

                % iterate through a scan and fill in L4 activity and bar's 
                % X/Y locations at each scan step
                for iScan = 1:stimParams.scanLength

                    % set coordinates of the bar and update xInit
                    % and yInit for the next step
                    if orientation == 1 % 0˚
                        xAtScanStep = xInit:1:(xInit + (stimParams.barLength - 1));
                        yAtScanStep = repelem(yInit, stimParams.barLength);

                        xInit = xInit;
                        yInit = yInit + direction;
                    elseif orientation == 2 % 45˚
                        xAtScanStep = xInit:1:(xInit + (stimParams.barLength - 1));
                        yAtScanStep = yInit:1:(yInit + (stimParams.barLength - 1));

                        % alternate between steps in x and y directions
                        xInit = xInit - (mod(iScan, 2) == 0) * direction;
                        yInit = yInit + (mod(iScan, 2) == 1) * direction;                      
                    elseif orientation == 3 % 90˚
                        xAtScanStep = repelem(xInit, stimParams.barLength);
                        yAtScanStep = yInit:-1:(yInit - (stimParams.barLength - 1));

                        xInit = xInit + direction;
                        yInit = yInit;
                    elseif orientation == 4 % 135˚
                        xAtScanStep = xInit:-1:(xInit - (stimParams.barLength - 1));
                        yAtScanStep = yInit:1:(yInit + (stimParams.barLength - 1));

                        % alternate between steps in x and y directions
                        xInit = xInit + (mod(iScan, 2) == 0) * direction;
                        yInit = yInit + (mod(iScan, 2) == 1) * direction;
                    end
    

                    if stimParams.isWraparound
                        % wrap around for bars that go over the edge of the
                        % receptive field
                        % set the x and y locations of the bar
                        % NOTE: this is fixed from IGOR code
                        wrapX = @(x) (1 + mod(x-1, stimParams.nX));
                        wrapY = @(y) (1 + mod(y-1, stimParams.nY));

                        xAtScanStep = wrapX(xAtScanStep);
                        yAtScanStep = wrapY(yAtScanStep);
                    else
                        % truncate the bar where it goes past the
                        % boundaries of the receptive field
                        invalidPoints = xAtScanStep <= 0 | yAtScanStep <= 0 | xAtScanStep > stimParams.nX | yAtScanStep > stimParams.nY;
                        xAtScanStep(invalidPoints) = NaN;
                        yAtScanStep(invalidPoints) = NaN;
                    end
                    barXLoc(:, iScan) = xAtScanStep;
                    barYLoc(:, iScan) = yAtScanStep;
    
                    
                    % TODO: documentation: describe L4activity array
                    % structure
                    
                    if stimParams.isFoldiak
                        % (Venkatesh L4 indexing)
                        % find index of the beginning of the set of L4 neurons corresponding to this bar's orientation
                        iOrientedL4 = (orientation - 1) * stimParams.nX * stimParams.nY;
                        % set the corresponding L4 neurons to ON
                        iL4 = iOrientedL4 + ((yAtScanStep - 1) * stimParams.nX + xAtScanStep);
                    else   
                        % (Berry L4 orientation indexing)
                        % NEED TO DEBUG
                        iL4 = orientation + stimParams.nOrient * ((yAtScanStep - 1) * stimParams.nX + (xAtScanStep - 1));
                    end


                    % DEBUG 2/17
                    iL4 = iL4(~isnan(iL4));
                    %
                    L4Activity(iL4, iScan) = 1;
                end
            end

            for iScan = 1:stimParams.scanLength
                if strcmp(stimParams.noiseType, 'exact')
                    numNoisyL4 = stimParams.noiseAmt;
                elseif strcmp(stimParams.noiseType, 'poisson')
                    numNoisyL4 = poissrng(stimParams.noiseAmt);
                end

                addedActiveL4 = randi(256, numNoisyL4, 1);
                while sum(L4Activity(addedActiveL4), 'all') > 0
                    addedActiveL4(L4Activity(addedActiveL4) == 1) = ...
                        randi(256, sum(L4Activity(addedActiveL4) == 1), 1);
                end          

                L4Activity(addedActiveL4, iScan) = 1;
            end

            % TODO: ADD THIS BACK IN
            % Add noise to scan.
            % 1. <numNoNoiseL4> neurons (out of 256) will NOT have uniform 
            %    noise added. For a non-pure-noise image, this includes the 
            %    bar plus a Poisson-distributed number of additional 'OFF' 
            %    L4 neurons.
            % 2. Of the noise-added neurons, set <nL4ExactNoise> to 'ON' (1).
            % numNoNoiseL4 = not(isNoiseImage) * stimParams.barLength + ...
            %     poissrnd(stimParams.nL4PoissNoise);
            % if numNoNoiseL4 > 256 numNoNoiseL4 = 256; end
            % 
            % if showOutput
            %     fprintf("Number of non-noisy L4 neurons = %d\n", numNoNoiseL4);
            % end
            % 
            % for iScan = 1:stimParams.scanLength
            % 
            %     [L4Sorted1, iSort1] = sort(L4Activity(:, iScan), 'descend'); % if a bar was created, put those entries at the front of the array
            % 
            %     iNoisyL4 = (1 + numNoNoiseL4):size(L4Sorted1, 1); % indices of noise-added L4 neurons
            %     L4Sorted1(iNoisyL4) = L4Sorted1(iNoisyL4) + 0.9 * rand(size(iNoisyL4)).'; % add noise to L4 neurons
            % 
            %     [L4Sorted2, iSort2] = sort(L4Sorted1, 'descend'); % resort to group the noisiest neurons on the right
            %     iSort2 = iSort1(iSort2); % track the indices of the L4 neurons
            % 
            %     iExactNoise = numNoNoiseL4:(numNoNoiseL4 + stimParams.nL4ExactNoise);
            %     L4Sorted2(iExactNoise) = 1; % set <nL4ExactNoise> of the L4 neurons to 1
            % 
            %     L4Activity(iSort2, iScan) = L4Sorted2; % unsort the L4 neurons to their original indices
            % end

            % SANJANA: add noise indiscriminately to neurons
            % L4Activity = L4Activity + 0 * randn(size(L4Activity)); % add Gaussian noise
            % L4Activity(L4Activity < 0) = 0;
            % L4Activity(L4Activity > 1) = 1;

            % TODO: return orientation?
        end
    end
end