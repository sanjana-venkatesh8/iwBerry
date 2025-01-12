classdef Foldiak
    %FOLDIAK Replica of Foldiak's model of a neuron with spatial invariance
    % http://dx.doi.org/10.1162/neco.1991.3.2.194
    %
    % Foldiak models a neuron that can identify objects from
    % different retinal position (i.e. as the retina drifts). These stimuli
    % are modeled as bars of length 4 that scan across an 8x8 grid (visual
    % field). 
    % NOTE: Foldiak uses bars of infinite length and scans in two 
    % directions. We use bars of length 4 that scan in 4 directions.
    
    properties
        stimuliSets StimuliParams
        connections;

        nComplex;
    end
    
    methods
        function [obj, initWeights] = Foldiak(args)
            %FOLDIAK Construct an instance of this class
            %   Detailed explanation goes here

            arguments
                args.stimuliSets StimuliParams
                
                args.weightsRange
                args.nComplex
            end

            obj.nComplex = args.nComplex;
            
            % Fully connect the simple neurons to each of the complex
            % units, initializing each weight randomly and uniformly on 
            % (0, 0.1).
            obj.stimuliSets = args.stimuliSets;

            nX = args.stimuliSets.nX;
            nY = args.stimuliSets.nY;
            nOrient = args.stimuliSets.nOrient;
            obj.connections = args.weightsRange * ...
                rand(nX * nY * nOrient, args.nComplex);
            initWeights = obj.connections;
        end
        
        function weights = train(obj, args)
            %TRAIN This function trains the four complex units.
            %   Bar stimuli are passed through the model and weights are
            %   updated according to Foldiak's modified Hebbian rule:
            %       ∆w_ij(t) = alpha * yTrace_i(t) * (x_j(t) - w_ij(t)) 

            arguments
                obj Foldiak
                args.delta
                args.alpha
            end
            
            stimuli = obj.stimuliSets.trainSet.L4Activity;
            activityTrace = zeros(1, obj.nComplex);

            for i = 1:size(stimuli, 3) % iterate through stimuli
                for j = 1:size(stimuli, 2) % iterate through each scan step in a stimulus
                    stimBar = stimuli(:, j, i);

                    % The complex unit with the highest weighted sum of inputs
                    % is activated.
                    weightedSum = sum(stimBar .* obj.connections, 1);
                    activity = (weightedSum == max(weightedSum));
    
                    % update average postsynaptic activity (trace) according to
                    % rule in Foldiak paper.
                    activityTrace = (1 - args.delta) * activityTrace + args.delta * activity; 
    
                    % update synaptic weights
                    % STOPPED HERE: ISSUE: trace not working. weight
                    % updates have a pattern

                    deltaWeight = args.alpha * activityTrace .* ...
                        (repmat(stimBar, 1, obj.nComplex) - obj.connections);
                    obj.connections = obj.connections + deltaWeight;
                    
                    weights = obj.connections;
                    % DEBUG
                    % pause(0.1)
                    % FoldiakGraphics.showConnections(8,8,4,4,obj.connections)
                    % h2 = subplot(2,1,2);
                    % clf(h2)
                    % barXLoc = obj.stimuliSets.trainSet.barXLoc(:, j, i);
                    % barYLoc = obj.stimuliSets.trainSet.barYLoc(:, j, i);
                    % % h1 = axes;
                    % scatter(barXLoc, barYLoc, "filled")
                    % set(h2,'YDir','reverse')
                    % xlim([1 8])
                    % ylim([1 8])
                    % 
                    % subplot(2,1,1);
                    % histogram(obj.connections);
                end
            end
        end
    end
end

