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
        nL4ExactNoise;
        nL4PoissNoise;
        propNoiseScan %double {mustBeInRange(propNoiseScan, 0, 1)}
        nTestRepeats %int32
        nTestInst %int32
        nTrainInst %int32
        nDendRecord %int32
        recordingTimeInterval;
        nL4Inputs;
        % non-input properties
        trainSet;
        testSet;
    end
    
    methods (Access = public)
        function stimParams = StimuliParams(nX, nY, ...
                nOrient, barLength, scanLength, nL4ExactNoise, ...
                nL4PoissNoise, propNoiseScan, nTestRepeats, nTestInst, ...
                nTrainInst, nDendRecord, recordingTimeInterval)
            
            %STIMULIPARAMS Construct an instance of this class

            % TODO: CHANGE TO NAME-VALUE ARGS
            % arguments
            % end

            stimParams.nX = nX;
            stimParams.nY = nY;
            
            stimParams.nOrient = nOrient;

            stimParams.barLength = barLength;

            stimParams.scanLength = scanLength;
            % mustBeMember(scanLength, 1:8) hard coded 8?

            stimParams.nL4ExactNoise = nL4ExactNoise;
            stimParams.nL4PoissNoise = nL4PoissNoise;
            stimParams.propNoiseScan = propNoiseScan;
            stimParams.nTestRepeats = nTestRepeats;
            stimParams.nTestInst = nTestInst;
            stimParams.nTrainInst = nTrainInst;
            stimParams.nDendRecord = nDendRecord;
            stimParams.recordingTimeInterval = recordingTimeInterval;
            % TODO: what to do if missing params in argument??

            % TODO: include showOutput param
            stimParams.nL4Inputs = nX * nY * nOrient;
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
            
            testSet.L4Activity = zeros(stimParams.nL4Inputs, ...
                stimParams.scanLength, stimParams.nTestInst);
            testSet.barXLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTestInst); 
            testSet.barYLoc = zeros(stimParams.barLength, ...
                stimParams.scanLength, stimParams.nTestInst); 

            % generate training set
            for iTrainInst = 1:stimParams.nTrainInst
                [trainSet.L4Activity(:, :, iTrainInst), ...
                    trainSet.barXLoc(:, :, iTrainInst), ...
                    trainSet.barYLoc(:, :, iTrainInst)] ...
                                        = makeScan(stimParams=stimParams, ...
                                                   showOutput=showOutput);
            end

            % generate testing set
            for iTestInst = 1:stimParams.nTestInst
                [testSet.L4Activity(:, :, iTestInst), ...
                    testSet.barXLoc(:, :, iTestInst), ...
                    testSet.barYLoc(:, :, iTestInst)] ...
                                        = makeScan(stimParams=stimParams, ...
                                                   showOutput=showOutput);
            end

        end

        % generate one scan and its corresponding L4 activity and X/Y
        % locations at each timestep
        function [L4Activity, barXLoc, barYLoc] = makeScan(args) % equivalent of Do_Bar_Scan
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

            isNoiseImage = rand() < stimParams.propNoiseScan;
            if not(isNoiseImage)
                direction = 2 * randi(2) - 3; % set random scan direction (-1 or 1)
                orientation = randi(stimParams.nOrient); % set random orientation of bar ([1, nOrient])

                if showOutput
                    fprintf("Orientation = %d, direction = %d\n", orientation, direction)
                end

                % set random initial bar location
                xInit = randi(stimParams.nX);
                yInit = randi(stimParams.nY);

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
    

                    % wrap around for bars that go over the edge of the
                    % receptive field
                    % set the x and y locations of the bar
                    % NOTE: this is fixed from IGOR code
                    wrapX = @(x) (1 + mod(x-1, stimParams.nX));
                    wrapY = @(y) (1 + mod(y-1, stimParams.nY));

                    xAtScanStep = wrapX(xAtScanStep);
                    yAtScanStep = wrapY(yAtScanStep);

                    barXLoc(:, iScan) = xAtScanStep;
                    barYLoc(:, iScan) = yAtScanStep;
    
                    
                    % TODO: documentation: describe L4activity array
                    % structure
                    
                    % (Venkatesh L4 indexing)
                    % find index of the beginning of the set of L4 neurons corresponding to this bar's orientation
                    % iOrientedL4 = (orientation - 1) * stimParams.nX * stimParams.nY;
                    % set the corresponding L4 neurons to ON
                    % iL4 = iOrientedL4 + ((yAtScanStep - 1) * stimParams.nX + xAtScanStep);

                    % (Berry L4 orientation indexing)
                    % NEED TO DEBUG
                    iL4 = orientation + stimParams.nOrient * ((yAtScanStep - 1) * stimParams.nX + (xAtScanStep - 1));

                    L4Activity(iL4, iScan) = 1;
                end
            end

            % TODO: ADD THIS BACK IN
            % Add noise to scan.
            % 1. <numNoNoiseL4> neurons will not have uniform noise added.
            %    For a non-pure-noise image, this includes the bar plus a
            %    Poisson-distributed number of additional 'OFF' L4 neurons.
            % 2. Of the noise-added neurons, set <nL4ExactNoise> to 'ON' (1).
            % numNoNoiseL4 = not(isNoiseImage) * stimParams.barLength + ...
            %     poissrnd(stimParams.nL4PoissNoise);
            % if showOutput
            %     fprintf("Number of non-noisy L4 neurons = %d\n", numNoNoiseL4);
            % end
            % 
            % for iScan = 1:stimParams.scanLength
            % 
            %     [L4Sorted1, iSort1] = sort(L4Activity(:, iScan), 'descend'); % if a bar was created, put those entries at the front of the array
            % 
            %     iNoisyL4 = (1 + numNoNoiseL4):size(L4Sorted1); % indices of noise-added L4 neurons
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

            % TODO: return orientation?
        end
    end
end