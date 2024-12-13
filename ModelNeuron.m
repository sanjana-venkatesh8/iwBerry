classdef ModelNeuron
%MODELNEURON Creates a model neuron and models plasticity in response to bar stimuli by passing it train and test stimuli.

    properties
        % TODO: VALIDATE?
        dendParams DendriteParams
        stimParams StimuliParams
        % non-input properties
        nL4Neurons;
        synL4Inputs;
        synUInput;
        synGInhib;
        synPlastWindow;
        nSynRecycles;
        synInputWMax;
        plasticityFlag
        areBranchesLeaky logical
    end
    
    % TODO: change access
    methods (Access = public)
        function modelNeuron = ModelNeuron(args)
        %ModelNeuron Initializes the model neuron. Part of Jens2_Plasticity.
            % TODO: fill in argument validation
            arguments
                args.dendParams DendriteParams
                args.stimParams StimuliParams
                args.plasticityFlag {mustBeMember(args.plasticityFlag, [1 2 3 4])}
            end

            modelNeuron.dendParams = args.dendParams;
            modelNeuron.stimParams = args.stimParams;
            modelNeuron.plasticityFlag = args.plasticityFlag;
            
            % Initialize the model neuron
            modelNeuron.nL4Neurons = modelNeuron.stimParams.nX * modelNeuron.stimParams.nY * modelNeuron.stimParams.nOrient;
            % randomly assign an L4 input neuron to each synapse
            modelNeuron.synL4Inputs = randi(modelNeuron.nL4Neurons, modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize);
            % randomly assign an input potential (u) in [-weightsRange, weightsRange] 
            % to each synapse
            modelNeuron.synUInput = modelNeuron.dendParams.weightsRange .* ...
                (2 .* rand(modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize) - 1); 
            % set initial inhibitory conductance for each synapse
            modelNeuron.synGInhib = modelNeuron.dendParams.initSynGInhib .* ...
                ones(modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize);
            % plasticity window for each synapse
            modelNeuron.synPlastWindow = zeros(modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize);
            % counter variable for number of recycling events per synapse
            modelNeuron.nSynRecycles = zeros(modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize);
            
            modelNeuron.areBranchesLeaky = (modelNeuron.dendParams.branchGLeak > 0);
            if ~modelNeuron.areBranchesLeaky
                % maximum EPSP (excitatory post-synaptic potential) for each synapse
                modelNeuron.synInputWMax = modelNeuron.dendParams.synStrength .* ...
                    ones(modelNeuron.dendParams.nBranches, modelNeuron.dendParams.branchSize);
            else
                % If branches are leaky, encode inhibitory conductance
                % TODO: understand why this is
                modelNeuron.synInputWMax = Constants.vExc ./ (1 + modelNeuron.dendParams.branchGLeak + modelNeuron.synGInhib);
            end
            
        end

        function [resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter, modelNeuron] ...
                = jens2Plasticity(modelNeuron)
            arguments
                modelNeuron ModelNeuron
            end

            fprintf("\n**************************************************\n" + ...
                      "***            STARTING SIMULATION             ***\n" + ...
                      "**************************************************\n\n");


            % CALCULATE STATISTICS BEFORE PLASTICITY, USING TEST SET
            % STIMULI
            [resultsBefore, modelNeuron] = stimulusUpdate(modelNeuron, plasticityFlag=1, isTrain=false);
            
            % TODO: where do 51 and 21 come from?
            UHistogram_Before = histcounts(modelNeuron.synUInput, NumBins=51);
            WMaxHistogram_Before = histcounts(modelNeuron.synInputWMax, NumBins=21);
            % figure("U Before")
            % histogram(modelNeuron.synUInput, NumBins=51);
            % figure("W_{max} Before")
            % histogram(modelNeuron.synInputWMax, NumBins=21);

            fprintf("Before plasticity: \n\ttotal # NMDA spikes = %d \n\ttotal # somatic spikes = %d\n", ...
                sum(resultsBefore.NMDASpikesPerTimestep), sum(resultsBefore.didSomaSpikePerTimestep));

            % analyze receptive fields
            resultsRFBefore = dendriticRFAnalyze(modelNeuron=modelNeuron, branchSpikeHist=resultsBefore.branchRF, ...
                cumulStimHist=resultsBefore.cumulStimRF, showOutput=true);

            nBranches = modelNeuron.dendParams.nBranches;
            resultsRFBefore.branchIOrient = resultsRFBefore.branchIOrient(1:nBranches) .* resultsBefore.branchNMDASpikeRate ...
                ./ sum(resultsBefore.branchNMDASpikeRate);
            totalIndex_Before = sum(resultsRFBefore.branchIOrient);
            resultsRFBefore.branchSize1 = resultsRFBefore.branchSize1(1:nBranches) .* resultsBefore.branchNMDASpikeRate ...
                ./ sum(resultsBefore.branchNMDASpikeRate);
            totalSize_Before = sum(resultsRFBefore.branchSize1);

            % IMPLEMENT SYNAPTIC PLASTICITY
            [resultsPlast, modelNeuron] = stimulusUpdate(modelNeuron, plasticityFlag=modelNeuron.plasticityFlag, isTrain=true);

            % SOMETHING WRONG HERE
            synRecycleHistogram = histcounts(modelNeuron.nSynRecycles, NumBins=int32(max(modelNeuron.nSynRecycles, [], 'all') + 1));
            % figure(Name="Number of Recycled Synapses")
            % histogram(modelNeuron.nSynRecycles, NumBins=int32(max(modelNeuron.nSynRecycles, [], 'all') + 1));

            tau = 30;
            x = 1:(10*tau) - 1;
            kernel = x .* exp(-x ./ tau);
            kernel = kernel / sum(kernel);

            convolvedDidSomaSpikePerTimestep = conv(kernel, resultsPlast.didSomaSpikePerTimestep);
            
            % CALCULATE STATISTICS AFTER PLASTICITY, ON TEST SET
            [resultsAfter, modelNeuron] = stimulusUpdate(modelNeuron, plasticityFlag=1, isTrain=false);

            nTimesteps = size(modelNeuron.stimParams.testSet.L4Activity, 3) * ...
                    size(modelNeuron.stimParams.testSet.L4Activity, 2);
            fprintf("After plasticity: \r\n\ttotal number of NMDA spikes = %d\r\n\ttotal number of somatic spikes = %d\r\n", ...
                sum(resultsAfter.NMDASpikesPerTimestep), sum(resultsAfter.didSomaSpikePerTimestep));
            
            % TODO: FIX THE WAY THIS IS PRINTING
            NMDA_SpikesPerTimestep_Rate = resultsAfter.NMDASpikesPerTimestep ./ nTimesteps;
            % fprintf("NMDA Spikes per timestep = %d\r\nNMDA spikes per branch = %d\r\n\n", ...
            %      NMDA_SpikesPerTimestep_Rate, NMDA_SpikesPerTimestep_Rate ./ modelNeuron.dendParams.nBranches);

            if (modelNeuron.dendParams.somaGLeak == 0)
                fprintf("After plasticity:\r\n\t soma threshold = %d\r\n\tmV spike rate = %d\r\n", ...
                    resultsPlast.somaThreshPerTimestep(end), sum(resultsAfter.didSomaSpikePerTimestep) / nTimesteps);
            else
                fprintf("After plasticity:\r\n\t soma gain = %d\r\n\tspike rate = %d\r\n", ...
                    resultsPlast.somaGainPerTimestep(end), sum(resultsAfter.didSomaSpikePerTimestep) / nTimesteps);
            end

            W_After = modelNeuron.synInputWMax ./ (1 + exp(-1 * modelNeuron.synUInput));

            % TODO: where do 51 and 21 come from?
            UHistogram_After = histcounts(modelNeuron.synUInput, NumBins=51);
            WMaxHistogram_After = histcounts(modelNeuron.synInputWMax, NumBins=21);
            % figure(Name="U After")
            % histogram(modelNeuron.synUInput, NumBins=51);
            % figure(Name="W_{max} after")
            % histogram(modelNeuron.synInputWMax, NumBins=21);

            % analyze receptive fields
            resultsRFAfter = dendriticRFAnalyze(modelNeuron=modelNeuron, branchSpikeHist=resultsAfter.branchRF, ...
                cumulStimHist=resultsBefore.cumulStimRF, showOutput=true);
            % TODO: ISSUE WITH CUMULSTIMRF_BEFORE

            resultsRFAfter.branchIOrient = resultsRFAfter.branchIOrient(1:nBranches) .* resultsAfter.branchNMDASpikeRate ...
                ./ sum(resultsAfter.branchNMDASpikeRate);
            totalIndex_After = sum(resultsRFAfter.branchIOrient);
            resultsRFAfter.branchSize1 = resultsRFAfter.branchSize1(1:nBranches) .* resultsAfter.branchNMDASpikeRate ...
                ./ sum(resultsAfter.branchNMDASpikeRate);
            totalSize_After = sum(resultsRFAfter.branchSize1);

            fprintf("Before plasticity: \r\n\tAverage orientation index = %d \r\n\tAverage size = %d\r\n", ...
                totalIndex_Before, totalSize_Before);
            fprintf("After plasticity: \r\n\tAverage orientation index = %d \r\n\tAverage size = %d\r\n", ...
                totalIndex_After, totalSize_After);
            

            fprintf("\n**************************************************\n" + ...
                      "***            SIMULATION FINISHED             ***\n" + ...
                      "**************************************************\n\n");
        end

        function [results, modelNeuron] = stimulusUpdate(modelNeuron, args)
        %stimulusUpdate Runs a set of stimuli through the model and calculates statistics.

        % equivalent of Do_Run_BarStim
        % From MJB: This function loops through a set of stimuli, calls 
        % 'Dendritic_Clusters' once per time step, and implements somatic
        % adaptation.

            arguments
                modelNeuron ModelNeuron
                args.plasticityFlag {mustBeMember(args.plasticityFlag, 1:4)}
                    % plasticityFlag == 1		no plasticity
                    % plasticityFlag == 2		somatic adaptation only
                    % plasticityFlag == 3		somatic adaptation + plasticity of dendritic branches
                    % plasticityFlag == 4		all plasticity (i.e. somatic + dendritic + synaptic)
                args.isTrain logical
            end

            % Calculate number of points to record at based on recording
            % interval and training set size.
            if args.isTrain
                barSetL4Activity = modelNeuron.stimParams.trainSet.L4Activity;
            else
                barSetL4Activity = modelNeuron.stimParams.testSet.L4Activity;
            end
            nBars = size(barSetL4Activity, 3);

            nTimesteps = nBars * ...
                    size(barSetL4Activity, 2);
            nRecordPoints = ceil(nTimesteps / modelNeuron.stimParams.recordingTimeInterval);
      
            if modelNeuron.stimParams.nDendRecord > 0
                synURecord = zeros(nRecordPoints, ...
                    modelNeuron.dendParams.branchSize, ...
                    modelNeuron.stimParams.nDendRecord);
                synL4Record = zeros(nRecordPoints, ...
                    modelNeuron.dendParams.branchSize, ...
                    modelNeuron.stimParams.nDendRecord);
                
                branchWMaxRecord = zeros(nRecordPoints, ...
                    modelNeuron.stimParams.nDendRecord);
                branchVRecord = zeros(nRecordPoints, ...
                    modelNeuron.stimParams.nDendRecord);
                branchSpikeRecord = zeros(nRecordPoints, ...
                    modelNeuron.stimParams.nDendRecord); % TODO: idk what this does
            end

            % TODO: how were these constants chosen? - need to assign
            % somabinsize somewhere
            numBins = 200;
            somaHistogram = zeros(2, numBins+1);
            if ~modelNeuron.areBranchesLeaky
                somaBinSize = 0.5;
            else
                somaBinSize = 0.2;
            end
            somaHistogram(1, :) = 0:somaBinSize:(somaBinSize*numBins);
        
            % Create tracker variables
            didSomaSpike = false;
            somaThresh = modelNeuron.dendParams.somaThresh;
            somaGain = 1;
            somaGInhib = modelNeuron.dendParams.initSomaGInhib;

            nTotalSpikes = 0; % count total number of spikes

            % track important values at each timestep
            NMDASpikesPerTimestep = zeros(nTimesteps, 1);
            somaActivityPerTimestep = zeros(nTimesteps, 1);
            somaEPSPsPerTimestep = zeros(nTimesteps, 1);
            didSomaSpikePerTimestep = zeros(nTimesteps, 1);
            somaThreshPerTimestep = zeros(nTimesteps, 1);
            somaGainPerTimestep = zeros(nTimesteps, 1);
            somaInhibPerTimestep = zeros(nTimesteps, 1);

            % track cumulative branch activity and spikes per synapse
            branchActivity = zeros(modelNeuron.dendParams.nBranches, 1);
            branchDidNMDASpike = zeros(modelNeuron.dendParams.nBranches, 1);

            % track receptive fields
            cumulStimRF = zeros(modelNeuron.stimParams.nX, ...
                         modelNeuron.stimParams.nY, modelNeuron.stimParams.nOrient);
            branchRF = zeros((modelNeuron.dendParams.nBranches + 1), modelNeuron.stimParams.nX, ...
                         modelNeuron.stimParams.nY, modelNeuron.stimParams.nOrient);

            % Iterate through all timesteps of the stimuli scans.

            iTime = 1;
            for iScan = 1:nBars
                for iStep = 1:modelNeuron.stimParams.scanLength
                    stimBarL4Activity = barSetL4Activity(:, iStep, iScan);
                    
                    % Calculate dendritic activity at this timestep.
                    [stepBranchActivity, stepBranchGExc, stepBranchGain, ...
                        stepBranchDidNMDASpike, stepTotalNMDASpikes, synActivity, modelNeuron] = ...
                        modelNeuron.calcDendActivity(...
                            doBackprop=didSomaSpike, ...
                            plasticityFlag=args.plasticityFlag, ...
                            stimBarL4Activity=stimBarL4Activity);

                     branchActivity = branchActivity + stepBranchActivity;
                     branchDidNMDASpike = branchDidNMDASpike + stepBranchDidNMDASpike;
                     % TODO: we might need uniform bin size across trials -
                     % i.e. change below line (bc it sets according to
                     % branchHist size)
                     % Also does it make sense to have this inside the loop?
                     branchHistogram = histcounts(stepBranchActivity, NumBins=ceil(modelNeuron.dendParams.branchThresh * 20));
                     % figure(Name="Branch activity")
                     % histogram(stepBranchActivity, NumBins=ceil(modelNeuron.dendParams.branchThresh * 20));

                     % Set gain at soma.
                     if ~modelNeuron.areBranchesLeaky
                         somaGain = 1;
                     else
                         somaGain = 1 / ...
                             (1 + modelNeuron.dendParams.somaGLeak + ...
                             somaGInhib);
                     end

                     % Calculate activity at soma due to L4 inputs.
                     somaNMDASpikes = sum(stepBranchDidNMDASpike) * ...
                         modelNeuron.dendParams.branchStrength * somaGain;
                     somaEPSPs = sum(sum(synActivity)) / ...
                         modelNeuron.dendParams.synAttenuation * somaGain;
                     somaActivity = somaNMDASpikes + somaEPSPs;

                     NMDASpikesPerTimestep(iTime) = sum(stepBranchDidNMDASpike);
                         % somaNMDASpikes / modelNeuron.dendParams.branchStrength / ...
                         % somaGain;
                     somaActivityPerTimestep(iTime) = somaActivity;
                     somaEPSPsPerTimestep(iTime) = somaEPSPs;
                     
                     didSomaSpike = (somaActivity >= somaThresh);
                     didSomaSpikePerTimestep(iTime) = didSomaSpike;
                     nTotalSpikes = nTotalSpikes + didSomaSpike; 

                     if args.plasticityFlag ~= 1
                         if ~modelNeuron.areBranchesLeaky
                             if didSomaSpike
                                 somaThresh = somaThresh + ...
                                     modelNeuron.dendParams.somaPlus;
                             else
                                 somaThresh = somaThresh + ...
                                     modelNeuron.dendParams.somaMinus;
                             end
                         else
                             if didSomaSpike
                                 % Increase inhibition to decrease somatic
                                 % gain.
                                 somaGInhib = somaGInhib + ...
                                     modelNeuron.dendParams.somaPlus;
                             else
                                 % Decrease inhibition to increase somatic
                                 % gain.
                                 somaGInhib = somaGInhib + ...
                                     modelNeuron.dendParams.somaMinus;
                             end
                             % Ensure somatic inhibition is nonnegative.
                             somaGInhib = max(somaGInhib, 0);
                         end
                     end % end plasticity flag

                     somaThreshPerTimestep(iTime) = somaThresh;
                     somaGainPerTimestep(iTime) = somaGain;
                     somaInhibPerTimestep(iTime) = somaGInhib;

                     % Equivalent of Dendritic_RF_Calc: Take the spike
                     % pattern in a single time bin and increment the
                     % histograms for cumulative stimulus and branch
                     % receptive field (RF).
                     
                     for iL4Input = 1:modelNeuron.nL4Neurons
                         % If this L4 neuron is active
                         if stimBarL4Activity(iL4Input) == 1 
                             % Not confident in this
                             % Find orientation, X and Y location in RF
                             % corresponding to this L4 neuron
                             nXYLoc = modelNeuron.stimParams.nX * modelNeuron.stimParams.nY;
                             
                             % Berry L4 indexing
                             orientation = mod(iL4Input, modelNeuron.stimParams.nOrient);
                             if orientation == 0
                                 orientation = modelNeuron.stimParams.nOrient;  % account for wraparound due to modulus
                             end

                             iX = mod((iL4Input - orientation) / modelNeuron.stimParams.nOrient, modelNeuron.stimParams.nX) ...
                                 + 1;
                             if iX == 0
                                 iX = modelNeuron.stimParams.nX; % account for wraparound due to modulus
                             end

                             iY = ((iL4Input - orientation) / modelNeuron.stimParams.nOrient - (iX - 1)) / ...
                                 modelNeuron.stimParams.nX;
                             if iY == 0
                                 iY = modelNeuron.stimParams.nY; % account for wraparound due to modulus
                             end

                             % Venkatesh L4 indexing
                             % orientation = ceil(iL4Input / nXYLoc); %
                             % 
                             % iY = ceil(mod(iL4Input, nXYLoc) / modelNeuron.stimParams.nX);
                             % 
                             % iX = mod(iL4Input, modelNeuron.stimParams.nX);
                             % if iX == 0
                             %     iX = modelNeuron.stimParams.nX; % account for wraparound due to modulus
                             % end

                             % increment cumulative stimulus histogram
                             cumulStimRF(iX, iY, orientation) = ...
                                 cumulStimRF(iX, iY, orientation) + 1;

                             for iBranch = 1:modelNeuron.dendParams.nBranches
                                 % Increment branch RF histogram if it had
                                 % an NMDA spike.
                                 if stepBranchDidNMDASpike(iBranch) == 1
                                     branchRF(iBranch, iX, iY, orientation) = ...
                                         branchRF(iBranch, iX, iY, orientation) + 1;
                                 end
                             end

                             if didSomaSpike
                                 % The extra row in branchRF encodes the
                                 % somatic RF. Increment it if there was a
                                 % somatic spike.
                                 branchRF(modelNeuron.dendParams.nBranches + 1, ...
                                     iX, iY, orientation) = ...
                                     branchRF(modelNeuron.dendParams.nBranches + 1, ...
                                     iX, iY, orientation) + 1;
                             end
                         end
                     end

                     % TODO: FIX. If at a recording point, record synapse info
                     if mod(iTime, modelNeuron.stimParams.recordingTimeInterval) == 0
                         if modelNeuron.stimParams.nDendRecord > 0
                             % TODO: what does branchWMaxRecord mean??
                             iRecordPoint = iTime / modelNeuron.stimParams.recordingTimeInterval;
                             synURecord(iRecordPoint, :, :) = modelNeuron.synUInput(1:modelNeuron.stimParams.nDendRecord, :)';
                             synL4Record(iRecordPoint, :, :) = modelNeuron.synL4Inputs(1:modelNeuron.stimParams.nDendRecord, :)';

                             branchWMaxRecord(iRecordPoint, :) = modelNeuron.synInputWMax(1:modelNeuron.stimParams.nDendRecord, 1)';
                             branchVRecord(iRecordPoint, :) = stepBranchActivity(1:modelNeuron.stimParams.nDendRecord)';
                             branchSpikeRecord(iRecordPoint, :) = stepBranchDidNMDASpike(1:modelNeuron.stimParams.nDendRecord)';
                         end
                     end

                     iTime = iTime + 1;
                end
            end

            branchActivity = branchActivity / iTime;

            somaHistogram(2, 1:(end-1)) = histcounts(somaActivityPerTimestep, somaHistogram(1, :));

            % TODO: why isn't this ever used?
            % Number of NMDA spikes and frequency of this many spikes
            % happening per timestep
            NMDASpikeMax = max(NMDASpikesPerTimestep);
            NMDASpikeHist = zeros(2, int32(NMDASpikeMax+2));
            NMDASpikeHist(1, :) = 0:1:(NMDASpikeMax+1); % set bin edges
            NMDASpikeHist(2, 1:(end-1)) = histcounts(NMDASpikesPerTimestep, NMDASpikeHist(1, :));

            if args.plasticityFlag ~= 1
                modelNeuron.dendParams.somaThresh = somaThresh;
                modelNeuron.dendParams.initSomaGInhib = somaGInhib;
            end
                     
            % TODO: FIX THIS
            branchNMDASpikeRate = branchDidNMDASpike / nTimesteps; % firing probability per time bin
            avgVoltageRespPerBranch = branchActivity / nTimesteps; % avg voltage response per branch

            % create an object to store results
            results = ModelResults(nTotalSpikes, NMDASpikesPerTimestep, ...
                somaActivityPerTimestep, somaEPSPsPerTimestep, didSomaSpikePerTimestep, ...
                branchNMDASpikeRate, avgVoltageRespPerBranch, ...
                branchHistogram, somaHistogram, branchRF, cumulStimRF, ...
                somaThreshPerTimestep, somaGainPerTimestep, somaInhibPerTimestep, ...
                synURecord, synL4Record, branchWMaxRecord, branchVRecord, branchSpikeRecord);                                
        end

        % TODO: figure out where to insert this statement. If neuron is not leaky, use direct changes on wMax for gain
        % control. Else, use inhibitory shunting-based gain control.
        function [allBranchActivity, allBranchGExc, allBranchGain, ...
                allBranchDidNMDASpike, nTotalNMDASpikes, synActivity, modelNeuron] = ...
                calcDendActivity(modelNeuron, args)
        %calcDendActivity Takes input from one bar stimulus and calculates the corresponding dendritic activity.   
        %   Equivalent to Dendritic_Clusters
        %   calcDendActivity returns:
        %       1. allBranchActivity: sum of synaptic potentials on each branch
        %       2. allBranchGExc: excitatory conductance of each branch
        %       3. allBranchGain: gain on each branch
        %       4. allBranchDidNMDASpike: tracks if a spike occurred on each branch
        %       5. nTotalNMDASpikes: overall total number of NMDA spikes
        %       6. synActivity: weighted L4 input at each synapse

            arguments
                modelNeuron ModelNeuron
                args.stimBarL4Activity % TODO: validate size
                args.doBackprop logical
                args.plasticityFlag {mustBeMember(args.plasticityFlag, 1:4)}
                    % plasticityFlag == 1		no plasticity
                    % plasticityFlag == 2		somatic adaptation only
                    % plasticityFlag == 3		somatic adaptation + plasticity of dendritic branches
                    % plasticityFlag == 4		all plasticity (i.e. somatic + dendritic + synaptic)
            end
            
            branchGNa = modelNeuron.dendParams.branchGLeak * modelNeuron.dendParams.backpropAP ...
                            / (Constants.vNa + modelNeuron.dendParams.backpropAP); % TODO: what is this?

            % Set synaptic weights based on whether branch is leaky.
            if ~modelNeuron.areBranchesLeaky
                % TODO: make local wmax var
                synWeight = modelNeuron.synInputWMax ./ ...
                    (1 + exp (-1 * modelNeuron.synUInput));
            else
                % Use inhibitory shunting-based gain control
                % TODO: find out what this means ^

                % decode the synaptic inhibitory conductance (was encoded in
                % synInputWMax)
                % TODO: make local wmax var?
                modelNeuron.synGInhib = ...
                    Constants.vExc ./ modelNeuron.synInputWMax ...
                    - 1 - modelNeuron.dendParams.branchGLeak;
                % synWeight encodes the conductance of each synapse
                synWeight = 1 ./ (1 + exp(-1 * modelNeuron.synUInput));
            end

            % % FOR BERRY TESTING ONLY
            % % a vertical bar with:
            % orientation = 3;
            % xAtScanStep = [4 4 4 4] + 1;
            % yAtScanStep = [1 2 3 4] + 1;
            % 
            % %Venkatesh L4 indexing
            % % % find index of the beginning of the set of L4 neurons corresponding to this bar's orientation
            % % iOrientedL4 = (orientation - 1) * modelNeuron.stimParams.nX * modelNeuron.stimParams.nY;
            % % % set the corresponding L4 neurons to ON
            % % iL4 = iOrientedL4 + ((yAtScanStep - 1) * modelNeuron.stimParams.nX + xAtScanStep);
            % 
            % iL4 = orientation + modelNeuron.stimParams.nOrient * ((yAtScanStep - 1) * modelNeuron.stimParams.nX + (xAtScanStep - 1));
            % 
            % args.stimBarL4Activity = zeros(modelNeuron.stimParams.nL4Inputs, 1);
            % args.stimBarL4Activity(iL4) = 1;
            % 
            % modelNeuron.synL4Inputs = readmatrix("berryDebug/Jens_input_L4_save0.txt");
            % modelNeuron.synL4Inputs = modelNeuron.synL4Inputs + 1; % correct off-by-one indexing
            % 
            % % modelNeuron.synL4Inputs = readmatrix("b2VL4Indices.txt");
            % modelNeuron.synUInput = readmatrix("berryDebug/Jens_input_u_save0.txt");
            % % synWeight = modelNeuron.synInputWMax ./ ...
            % %         (1 + exp (-1 * modelNeuron.synUInput)); % TODO: THIS IS FOR NONLEAKY
            % 
            % modelNeuron.synGInhib = ...
            %         Constants.vExc ./ modelNeuron.synInputWMax ...
            %         - 1 - modelNeuron.dendParams.branchGLeak;
            % % synWeight encodes the conductance of each synapse
            % synWeight = 1 ./ (1 + exp(-1 * modelNeuron.synUInput));
            % % END BERRY TESTING BLOCK

            % Set synapse activity based on which L4 neurons are turned ON
            % by the input.
            synActivity = args.stimBarL4Activity(modelNeuron.synL4Inputs);
            % Multiply ON/OFF input by synaptic weight.
            synActivity = synActivity .* synWeight;
            % Add multiplicative white Gaussian noise to each synapse
            synActivity = synActivity .* ...
                (1 + modelNeuron.dendParams.synNoise * randn(modelNeuron.dendParams.nBranches, ...
                modelNeuron.dendParams.branchSize));

            % Create output variables
            allBranchActivity = zeros(modelNeuron.dendParams.nBranches, 1);
            allBranchGExc = zeros(modelNeuron.dendParams.nBranches, 1);
            allBranchGain = zeros(modelNeuron.dendParams.nBranches, 1);
            allBranchDidNMDASpike = zeros(modelNeuron.dendParams.nBranches, 1);
            nTotalNMDASpikes = 0;

            for iBranch = 1:modelNeuron.dendParams.nBranches
                % sum activity along all synapses in the branch
                branchGExc = sum(synActivity(iBranch, :));

                branchSynGInhib = modelNeuron.synGInhib(iBranch, 1);

                if ~modelNeuron.areBranchesLeaky
                    % Use direct changes on wMax for gain control
                    % If doBackprop is true, add de
                    % polarization due to a
                    % backpropagating action potential
                    branchActivity = branchGExc + ...
                        (args.doBackprop == true) * modelNeuron.dendParams.backpropAP;
                    branchGain = 1;
                else
                    % Use inhibitory shunting-based gain control.
                    branchActivity = branchGExc * Constants.vExc + ...
                        (args.doBackprop == true) * branchGNa * Constants.vNa;
                   
                    denom = branchGExc + modelNeuron.dendParams.branchGLeak + ...
                        branchSynGInhib + (args.doBackprop == true) * branchGNa;
                    branchActivity = branchActivity ./ denom;
                    branchGain = Constants.vExc ./ denom;
                end

                % Apply inhibitory gain control to synapse activity (EPSPs).
                synActivity(iBranch,:) = synActivity(iBranch,:) .* branchGain;

                % Calculate actual spike threshold and determine if a spike
                % occurred.
                actualThresh = modelNeuron.dendParams.branchThresh * ...
                    (1 + modelNeuron.dendParams.branchNoise * randn);
                didBranchSpike = branchActivity >= actualThresh;

                % Update storage variables
                allBranchActivity(iBranch)        = branchActivity;
                allBranchGExc(iBranch)            = branchGExc;
                allBranchGain(iBranch)            = branchGain;
                allBranchDidNMDASpike(iBranch)    = didBranchSpike;

                if args.plasticityFlag >= 3 % somatic + dendritic plasticity
                    % Adjust dendritic gain to change NMDA spike rate.
                    if didBranchSpike
                        if ~modelNeuron.areBranchesLeaky
                            % Scale synapse wMax to adjust NMDA spike
                            % rate
                            % TODO: make local wmax var?
                            modelNeuron.synInputWMax(iBranch, :) = ...
                                modelNeuron.synInputWMax(iBranch, :) .* ...
                                (1 + modelNeuron.dendParams.scaleNMDA);
                        else
                            % *Increment* (scaleNMDA < 0) inhibitory
                            % conductance to decrease gain.
                            modelNeuron.synGInhib(iBranch, :) = ...
                                modelNeuron.synGInhib(iBranch, :) - ...
                                modelNeuron.dendParams.scaleNMDA;
                        end
                    else
                        if ~modelNeuron.areBranchesLeaky
                            % TODO: make local wmax var?
                            modelNeuron.synInputWMax(iBranch, :) = ...
                                modelNeuron.synInputWMax(iBranch, :) .* ...
                                (1 + modelNeuron.dendParams.scaleNoNMDA);
                        else
                            % Decrement inhibitory conductance
                            % (scaleNoNMDA > 0)
                            modelNeuron.synGInhib(iBranch, :) = ...
                                modelNeuron.synGInhib(iBranch, :) - ...
                                modelNeuron.dendParams.scaleNoNMDA;
                        end
                    end
                end
                if args.plasticityFlag == 4 % add in synaptic plasticity
                    if didBranchSpike
                        % Set up plasticity time window.
                        modelNeuron.synPlastWindow(iBranch, :) = ...
                            modelNeuron.dendParams.plastTime;
                    else
                        modelNeuron.synPlastWindow(iBranch, :) = ...
                            modelNeuron.synPlastWindow(iBranch, :) - 1;
                    end

                    for iSyn = 1:modelNeuron.dendParams.branchSize
                        isSynActive = (synActivity(iBranch, iSyn) == 1);
                        if modelNeuron.synPlastWindow(iBranch, iSyn) > 0
                            modelNeuron.synUInput(iBranch, iSyn) = ...
                                modelNeuron.synUInput(iBranch, iSyn) + ...
                                ~isSynActive * modelNeuron.dendParams.duPotent + ...
                                isSynActive * modelNeuron.dendParams.duDepress;
                        else
                            modelNeuron.synUInput(iBranch, iSyn) = ...
                                modelNeuron.synUInput(iBranch, iSyn) + ...
                                ~isSynActive * modelNeuron.dendParams.duDecay + ...
                                isSynActive * modelNeuron.dendParams.duBaseline;
                        end

                        % Recycle spine (retract it and grow a new one)
                        % if synaptic potential is too low
                        if (modelNeuron.synUInput(iBranch, iSyn) < ...
                                modelNeuron.dendParams.uRecycle)

                            % Randomly assign a new input L4 neuron
                            modelNeuron.synL4Inputs(iBranch, iSyn) = ...
                                randi(256);
                            % Randomly set a new input potential
                            modelNeuron.synUInput(iBranch, iSyn) = ...
                                modelNeuron.dendParams.weightsRange .* ...
                                (2 .* rand() - 1); 
                            % Update the count of recycled synapses
                            modelNeuron.nSynRecycles(iBranch, iSyn) = ...
                                modelNeuron.nSynRecycles(iBranch, iSyn) + 1;
                        end
                    end
                end

                if didBranchSpike
                    nTotalNMDASpikes = nTotalNMDASpikes + 1;
                    % Reset ("shunt"?) the linear EPSPs when there is an
                    % NMDA spike.
                    synActivity(iBranch, :) = zeros(1, modelNeuron.dendParams.branchSize);
                end

            end % end iteration through branches

            % Cap synapse potentials at uMax.
            modelNeuron.synUInput = min(modelNeuron.synUInput, modelNeuron.dendParams.uMax);
            % Make sure there are no negative inhibitory conductances.
            modelNeuron.synGInhib = max(modelNeuron.synGInhib, 0);

            % Re-encode inhibitory conductance into synInputWMax
            if modelNeuron.areBranchesLeaky
                % TODO: make local wmax var?
                modelNeuron.synInputWMax = Constants.vExc ./ ...
                    (1 + modelNeuron.dendParams.branchGLeak + modelNeuron.synGInhib);
            end
        end
        
        function resultsRF = dendriticRFAnalyze(args) 
            %dendriticRFAnalyze Takes the histogram of stimuli causing an NMDA spike for each branch and calculates:
            % 1) orientation tuning curves for each branch
            % 2) spatial receptive fields for each branch
            % 'branch_spike_hist' is a 4D tensor with counts for [branch][x][y][orientation]
            % row 'num_branches + 1' is for the RF of the soma
            % creates output waves: 	'temp_branch_spatial_all'
            % 						    'temp_branch_orient_all'
            % 						    'temp_branch_size1'  &  'temp_branch_size2'
            % 						    'temp_branch_index'
            % 					    	'temp_branch_pref'
            % 					    	'temp_branch_vector'
            % Note: last two rows are: i) soma RF, ii) composite RF

            arguments
                args.modelNeuron ModelNeuron
                args.showOutput logical
                args.branchSpikeHist
                args.cumulStimHist
            end
    

           % assign abbreviated var names
           nBranches = args.modelNeuron.dendParams.nBranches;
           nX = args.modelNeuron.stimParams.nX;
           nY = args.modelNeuron.stimParams.nY;
           nL4Neurons = args.modelNeuron.nL4Neurons;
           nOrient = args.modelNeuron.stimParams.nOrient;


            if args.showOutput
                % fprintf("\nNumber of branches = %d\n" + ...
                %     "Number of x locations = %d\n", nBranches, nX);
                fprintf("\nNumber of branches = %d\n" + ...
                    "Number of x locations = %d\n" + ...
                    "Number of y locations = %d\n" + ...
                    "Number of orientations = %d\n" + ...
                    "Number of L4 neurons = %d\n", ...
                    nBranches, nX, nY, nOrient, nL4Neurons);
            end

            % sort the total count of stimuli by orientation and spatial location
            cumulOrient = zeros(nOrient, 1);
            cumulSpatial = zeros(nX, nY);

            for iOrient = 1:nOrient
                for iX = 1:nX
                    for iY = 1:nY
                        cumulOrient(iOrient) = cumulOrient(iOrient) + ...
                            args.cumulStimHist(iX, iY, iOrient);
                        cumulSpatial(iX, iY) = cumulSpatial(iX, iY) + ...
                            args.cumulStimHist(iX, iY, iOrient);
                    end
                end
            end
            
            compositeOrient = zeros(nOrient, 1);
            compositeSpatial = zeros(nX, nY);
            allBranchOrient = zeros(nBranches + 1, nOrient);
            allBranchSpatial = zeros(nBranches + 1, nX, nY);
            branchRFFrac = zeros(nOrient, nX, nY);
            branchRMax = zeros(nBranches + 1, 1);
            branchSize1 = zeros(nBranches + 1, 1);
            branchSize2 = zeros(nBranches + 1, 1);
            branchIOrient = zeros(nBranches + 1, 1);
            branchPref = zeros(nBranches + 1, 1);
            branchVector = zeros(nBranches + 1, 2);

            for iBranch = 1:nBranches
                orientTuning = zeros(nOrient, 1);
                spatialRF = zeros(nX, nY);

                % compile orientation tuning curve over all positions and
                % orientations
                
                for iOrient = 1:nOrient
                    orientTuning(iOrient) = sum(args.branchSpikeHist(iBranch,:,:,iOrient), 'all');
                end
                for iX = 1:nX
                    for iY = 1:nY 
                        spatialRF(iX, iY) = sum(args.branchSpikeHist(iBranch, iX, iY, :), 'all');
                    end
                end
                allBranchOrient(iBranch, :) = orientTuning ./ cumulOrient;
                allBranchSpatial(iBranch, :, :) = spatialRF ./ cumulSpatial;

                % only include dendritic branches in composite RF (last
                % row is soma?
                compositeOrient = compositeOrient + orientTuning;
                compositeSpatial = compositeSpatial + spatialRF;
                % if iBranch < nBranches
                %     compositeOrient = compositeOrient + orientTuning;
                %     compositeSpatial = compositeSpatial + spatialRF;
                % end

                branchRFFrac = args.branchSpikeHist(iBranch, :, :, :) ./ ...
                    reshape(args.cumulStimHist, size(args.branchSpikeHist(iBranch, :, :, :)));

                branchRMax(iBranch) = max(branchRFFrac,[],"all");
            end

            compositeOrient = compositeOrient ./ (args.modelNeuron.dendParams.nBranches - 1);
            compositeSpatial = compositeSpatial ./ (args.modelNeuron.dendParams.nBranches - 1);

            allBranchOrient(nBranches + 1, :) = ...
                compositeOrient ./ cumulOrient;
            allBranchSpatial(nBranches + 1, :, :) = ...
                compositeSpatial ./ cumulSpatial;

            % % DEBUG/CHECK - remove nan values
            % allBranchOrient(isnan(allBranchOrient)) = 0;
            % allBranchSpatial(isnan(allBranchSpatial)) = 0;

            for iBranch = 1:(nBranches + 1)
                spatialRF = allBranchSpatial(iBranch, :, :);
                
                RFSize1 = sum(spatialRF, 'all') / max(spatialRF, [], 'all');
                % if isnan(RFSize1) RFSize1 = 0; end % TODO: should this ever be NaN
                branchSize1(iBranch) = RFSize1;

                spatialRF = spatialRF ./ max(spatialRF, [], 'all');
                spatialRF = spatialRF > 0.2;
                RFSize2 = sum(sum(spatialRF));
                % if isnan(RFSize2) RFSize2 = 0; end % TODO: should this ever be NaN
                branchSize2(iBranch) = RFSize2;

                orientTuning = allBranchOrient(iBranch, :);
                pref = find(orientTuning == max(orientTuning), 1, 'first');
                null = find(orientTuning == min(orientTuning), 1, 'first');

                orientIndex = (orientTuning(pref) - orientTuning(null)) / ...
                    sum(orientTuning);
                % if isnan(orientIndex) orientIndex = 0; end % TODO: should this ever be NaN
                branchIOrient(iBranch) = orientIndex;
                branchPref(iBranch) = pref;
                
                orientNorm = rms(orientTuning);
                orientVecX = (orientTuning(2) - orientTuning(4)) / orientNorm;
                orientVecY = (orientTuning(1) - orientTuning(3)) / orientNorm;
                % if isnan(orientVecX) orientVecX = 0; end % TODO: should this ever be NaN
                % if isnan(orientVecY) orientVecY = 0; end % TODO: should this ever be NaN
                branchVector(iBranch, 1) = orientVecX;
                branchVector(iBranch, 2) = orientVecY;
            end

            % create an object to store results
            resultsRF = RFResults(allBranchSpatial, ...
                allBranchOrient, branchIOrient, branchPref, ...
                branchVector, branchSize1, branchSize2);
        end
    end
end