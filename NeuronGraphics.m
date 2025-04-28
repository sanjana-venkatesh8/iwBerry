classdef NeuronGraphics
    properties
        modelNeuronObj ModelNeuron

        resultsBefore ModelResults
        resultsPlast ModelResults
        resultsAfter ModelResults
        
        resultsRFBefore RFResults
        resultsRFAfter RFResults
    end

    methods
        % TODO: make these static

        function neuronGraphics = NeuronGraphics(...
                modelNeuronObj, resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter)
            arguments
                modelNeuronObj ModelNeuron

                resultsBefore ModelResults
                resultsPlast ModelResults
                resultsAfter ModelResults
                
                resultsRFBefore RFResults
                resultsRFAfter RFResults
            end

            neuronGraphics.modelNeuronObj = modelNeuronObj;

            neuronGraphics.resultsBefore = resultsBefore;
            neuronGraphics.resultsPlast = resultsPlast;
            neuronGraphics.resultsAfter = resultsAfter;
            
            neuronGraphics.resultsRFBefore = resultsRFBefore;
            neuronGraphics.resultsRFAfter = resultsRFAfter;
        end

        % TODO: make these static methods ??
        function layout1(nG)
        %LAYOUT1 Plots orientation tuning vs. NMDA spike rate, colored by 
        % RF size, before and after plasticity.
            % TODO: might need * 100 for branchIOrient
            figure(Name="Layout 1: Orientation Tuning Across Subunits")            
            subplot(2, 1, 1)
            scatter(nG.resultsBefore.branchNMDASpikeRate, ...
                nG.resultsRFBefore.branchIOrient(1:(end-2)), 50, nG.resultsRFBefore.branchSize1(1:(end-2)), 'filled');
            c = colorbar;
            clim([0 15])
            c.Label.String = "RF Size (pixels)";
            title("Before plasticity")
            ylabel("Orientation Tuning Index")
            xlabel("NMDA Spike Rate")

            subplot(2, 1, 2)
            scatter(nG.resultsAfter.branchNMDASpikeRate, ...
                nG.resultsRFAfter.branchIOrient(1:(end-2)), 50, nG.resultsRFAfter.branchSize1(1:(end-2)), 'filled');
            c = colorbar;
            clim([0 15])
            c.Label.String = "RF Size (pixels)";
            title("After plasticity")
            ylabel("Orientation Tuning Index")
            xlabel("NMDA Spike Rate")
        end
        
        function layout2(nG, iBranch)
        %layout2 Plots the receptive fields of a given branch
            figure(Name="Layout 2: Branch Receptive Fields")
            t = tiledlayout('flow');

            nexttile
            % TODO: where to find EPSP and response?
            infoText = sprintf("Dendritic branch #%d\n" + ...
                                "Orientation index = %1.2f\n" + ...
                                "Receptive field size = %1.2f\n" + ...
                                "Max EPSP = %1.2f mV\n" + ...
                                "Max response = %d", ...
                                iBranch, nG.resultsRFAfter.branchIOrient(iBranch), ...
                                nG.resultsRFAfter.branchSize1(iBranch), ...
                                nG.modelNeuronObj.synInputWMax(iBranch, 1), NaN);
            text(0,0.7, infoText, FontSize=14)
            Ax = gca;
            Ax.Visible = 0;

            % plot branch receptive field
            spatialRF = squeeze(nG.resultsRFAfter.allBranchSpatial(iBranch, :, :)).';
            nexttile
            h = heatmap(spatialRF, 'CellLabelColor','none');
            h.Title = sprintf("Receptive field of Branch #%d after plasticity", iBranch);
            h.XLabel = "X location";
            h.YLabel = "Y location";
            colormap('hot')

            nexttile([1 2])
            nOrient = size(nG.resultsRFBefore.allBranchOrient, 2);
            nBranchSpikePerOrient = zeros(nOrient, 1);
            for i = 1:nOrient
                nBranchSpikePerOrient(i) = sum(nG.resultsAfter.branchRF(iBranch, :, :, i), 'all');
            end
            bar(nBranchSpikePerOrient)
            % histogram(nBranchSpikePerOrient, BinEdges=0.5+(0:1:nOrient))
            xticks(0:1:nOrient)
            title("Histogram of orientation of stimuli that caused a spike")
            xlabel("Bar orientation")
            ylabel("Count")
        end

        function layout3(nG)
        %LAYOUT3 Plots the soma receptive field vs. composite receptive
        %field
            % NEED TO FIX
            figure(Name="Layout 3: Soma Receptive Field vs. Composite Receptive Field")
            set(gcf, 'Position', [0 0 450 900])
            t = tiledlayout('flow');

            % Plot details for soma
            iSoma = size(nG.resultsRFBefore.allBranchSpatial, 1) - 1;
            nexttile
            % TODO: where to find EPSP and response?
            infoText = sprintf("Somatic receptive field\n" + ...
                                "Orientation index = %1.2f\n" + ...
                                "Receptive field size = %1.2f\n" + ...
                                "Max response = %d", ...
                                nG.resultsRFAfter.branchIOrient(iSoma), ...
                                nG.resultsRFAfter.branchSize1(iSoma), NaN);
            text(0,0.7, infoText, FontSize=14)
            Ax = gca;
            Ax.Visible = 0;

            % plot branch receptive field
            spatialRF = squeeze(nG.resultsRFAfter.allBranchSpatial(iSoma, :, :)).';
            nexttile
            h = heatmap(spatialRF, 'CellLabelColor','none');
            h.Title = "Receptive field of soma after plasticity";
            h.XLabel = "X location";
            h.YLabel = "Y location";
            colormap('hot')

            nexttile([1 2])
            nOrient = size(nG.resultsRFBefore.allBranchOrient, 2);
            nBranchSpikePerOrient = zeros(nOrient, 1);
            for i = 1:nOrient
                nBranchSpikePerOrient(i) = sum(nG.resultsAfter.branchRF(iSoma, :, :, i), 'all');
            end
            bar(nBranchSpikePerOrient)
            xticks(0:1:nOrient)
            title("Histogram of orientation of stimuli that caused a somatic spike")
            xlabel("Bar orientation")
            ylabel("Count")

            % Plot details for composite RF
            iComp = iSoma + 1;
            nexttile
            % TODO: where to find EPSP and response?
            infoText = sprintf("Composite receptive field\n" + ...
                                "Orientation index = %1.2f\n" + ...
                                "Receptive field size = %1.2f\n" + ...
                                "Max response = %d", ...
                                nG.resultsRFAfter.branchIOrient(iComp), ...
                                nG.resultsRFAfter.branchSize1(iComp), NaN);
            text(0,0.7, infoText, FontSize=14)
            Ax = gca;
            Ax.Visible = 0;

            % plot branch receptive field
            spatialRF = squeeze(nG.resultsRFAfter.allBranchSpatial(iComp, :, :)).';
            nexttile
            h = heatmap(spatialRF, 'CellLabelColor','none');
            h.Title = "Composite receptive field after plasticity";
            h.XLabel = "X location";
            h.YLabel = "Y location";
            colormap('hot')

            % TODO: make this step (and all histograms above) part of
            % DendriticRFAnalyze()
            nexttile([1 2])
            nOrient = size(nG.resultsRFBefore.allBranchOrient, 2);
            nBranchSpikePerOrient = zeros(nOrient, 1);
            for i = 1:nOrient
                nBranchSpikePerOrient(i) = sum(nG.resultsAfter.branchRF(1:(iSoma - 1), :, :, i), 'all');
            end
            bar(nBranchSpikePerOrient)
            xticks(0:1:nOrient)
            title("Composite histogram of orientation of stimuli that caused a spike")
            xlabel("Bar orientation")
            ylabel("Count")
        end
        
        function layout4(nG, iBranch, branchGLeak, tau)
        %LAYOUT4 Plots neuron adaptation: NMDA rate, synaptic max weight, and branch
        %inhibition over time

        % tau = smoothing window
        % TODO: rateTarget = NMDA firing rate target
        % TODO: fix x-lim of top graph

        % Calculate NMDA spike rate
        % Smooth spike train into a firing rate
        NMDASpikeRate = nG.resultsPlast.branchSpikeRecord(:, iBranch);
        % make a convolutional kernel with width tau
        x = 1:(10*tau) - 1;
        kernel = x .* exp(-x ./ tau);
        kernel = kernel / sum(kernel);
        NMDASpikeRate = conv(kernel, NMDASpikeRate);

        % Calculate branch inhibition
        % TODO: where does 60 come from?
        branchInhibPerTimestep = 60 ./ nG.resultsPlast.branchWMaxRecord - 1 - branchGLeak;

        title4 = sprintf("Layout 4: Adaptation of NMDA rate, Max weight, and branch inhibition in Branch #%d", ...
            iBranch);
        f = figure(Name=title4);
        subplot(3,1,1)
        % TODO: this is wrong
        plot(NMDASpikeRate);
        title(title4)
        ylabel("NMDA rate (Hz)")
        subplot(3,1,2)
        plot(nG.resultsPlast.branchWMaxRecord(:, iBranch));
        ylabel("Max weight (mV)")
        subplot(3,1,3)
        plot(branchInhibPerTimestep(:, iBranch));
        ylabel("Branch inhibition")
        xlabel("Recording step")
        end

        function layout5(nG, iBranch, iSyn, branchGLeak)
        %LAYOUT5 Plots internal potential, synapse's L4 input, and synaptic
        %weight
            % TODO: figure out how to calculate location based on synapse
            % ID
            L4FinalInput = nG.resultsPlast.synL4Record(end, iSyn, iBranch);
            orientation = mod(L4FinalInput, nG.modelNeuronObj.stimParams.nOrient);
             if orientation == 0
                 orientation = nG.modelNeuronObj.stimParams.nOrient;  % account for wraparound due to modulus
             end
    
             iX = mod((L4FinalInput - orientation) / nG.modelNeuronObj.stimParams.nOrient, nG.modelNeuronObj.stimParams.nX) ...
                 + 1;
             if iX == 0
                 iX = nG.modelNeuronObj.stimParams.nX; % account for wraparound due to modulus
             end
    
             iY = ((L4FinalInput - orientation) / nG.modelNeuronObj.stimParams.nOrient - (iX - 1)) / ...
                 nG.modelNeuronObj.stimParams.nX;
             if iY == 0
                 iY = nG.modelNeuronObj.stimParams.nY; % account for wraparound due to modulus
             end

            % Calculate synaptic weights
            branchExc = 1 ./ (1 + exp(-1 .* nG.resultsPlast.synURecord(:, iSyn, iBranch)));
            % TODO: branchInh might be wrong
            % TODO: where is 60 from?
            branchInh = 60 ./ nG.resultsPlast.branchWMaxRecord(iBranch) - 1 - branchGLeak;
            synWeights = 60 * branchExc / (branchExc + branchInh + branchGLeak);

            title5 = sprintf(...
                "Layout 5: Plasticity in Synapse %d on Branch %d\n" ...
                + "(final L4 input: x = %d, y = %d, orientation = %d)", ...
                iSyn, iBranch, iX, iY, orientation);
            figure(Name=title5)

            subplot(3,1,1)
            plot(synWeights)
            ylabel("Synaptic weight (mV)")
            title(title5)
            subplot(3,1,2)
            plot(nG.resultsPlast.synL4Record(:, iSyn, iBranch))
            ylabel("Input identity")
            subplot(3,1,3)
            plot(nG.resultsPlast.synURecord(:, iSyn, iBranch))
            ylabel("Internal potential")
            xlabel("Recording step")     
        end
        
        function layout7(nG)
            % TODO: actually write this
            figure(Name="Histograms over internal potential, maximum synapse strength, synaptic recycles")
            subplot(3,1,1)
            hold on
            [uBeforeCounts, uBeforeEdges] = histcounts(nG.resultsBefore.synURecord(1,:,:), NumBins=51);
            plot(uBeforeEdges(1:(end-1)), uBeforeCounts)
            [uAfterCounts, uAfterEdges] = histcounts(nG.resultsAfter.synURecord(1,:,:), NumBins=51);
            plot(uAfterEdges(1:(end-1)), uAfterCounts)
            % histogram(nG.resultsBefore.synURecord)
            % histogram(nG.resultsAfter.synURecord)
            title(sprintf("Histogram over values of internal potential (u) for %d branches", ...
                size(nG.resultsBefore.synURecord, 3)))
            xlabel("Internal potential")
            ylabel("Counts")
            legend("Before", "After")

            subplot(3,1,2)
            histogram(nG.modelNeuronObj.synInputWMax) % TODO: should record this in results
            title("Histogram over values of max possible synapse strength, w_{max}")
            xlabel("Value of w_{max}")
            ylabel("Number of synapses")

            subplot(3,1,3)
            histogram(nG.modelNeuronObj.nSynRecycles)
            title("Histogram over number of recycles per synapse")
            xlabel("Number of synaptic recycles")
            ylabel("Counts")
        end

        function plotStimulus(nG, args)
            arguments
                nG NeuronGraphics
                args.isTrain logical
                args.iStim
                args.isPlotWithOrientation
            end

            stimParams = nG.modelNeuronObj.stimParams;

            if args.isTrain
                stimSet = stimParams.trainSet;
                figure(Name=sprintf("Stimulus #%d in training set", args.iStim))
            else
                stimSet = stimParams.testSet;
                figure(Name=sprintf("Stimulus #%d in test set", args.iStim))
            end

            nTotStim = size(stimSet.L4Activity, 3);
            if args.iStim >= nTotStim || args.iStim < 1
                throw(MException('NeuronGraphics:IndexOutofBounds', ...
                    'Stimulus index out of bounds'))
            end

            L4Activity = stimSet.L4Activity(:, :, args.iStim);
            barXLoc = stimSet.barXLoc(:, :, args.iStim);
            barYLoc = stimSet.barYLoc(:, :, args.iStim);


            for iScan = 1:stimParams.scanLength
                if args.isPlotWithOrientation
                    for iOrient = 1:stimParams.nOrient
                        L4ChooseOrient = L4Activity(iOrient:4:256, iScan);
                        L4Heat = reshape(L4ChooseOrient, 8, 8).';
                        subplot(2, 2, iOrient)
                        hm = heatmap(L4Heat());
                        set(hm, 'ColorLimits', [0 1])
                        hm.Title = sprintf("L4 neurons with orientation %d", iOrient);
                        % colormap('hot')
                    end
                else
                    % another way of plotting the stimuli
                    h1 = axes;                        
                        scatter(barXLoc(:,iScan), barYLoc(:,iScan), "filled")
                    set(h1,'YDir','reverse')
                    xlim([1 8])
                    ylim([1 8])
                end
                pause(0.5)
            end
        end
        
        function branchSpikeRasterPlots(nG)
            %BRANCHSPIKERASTERPLOTS plots spiking patterns for branches and
            %soma. Adapted from Felix Schneider
            %(https://www.youtube.com/watch?v=27Y2c596-U0)
            figure(Units='normalized', Position=[0 0 0.3 1]); hold on
            subplot(2,1,1); hold on
            title("Spike firing for each branch")
            ylabel("Branch number")
            xlabel("Timestep")

            for iBranch = 1:nG.modelNeuronObj.dendParams.nBranches
                spikeTimes = find(nG.resultsAfter.branchSpikeRecord(:, iBranch).' == 1);
                xspikes = repmat(spikeTimes, 3, 1);
                yspikes = nan(size(xspikes));
             
                if ~isempty(yspikes)
                    yspikes(1, :) = iBranch - 1;
                    yspikes(2, :) = iBranch;
                end
             
                plot(xspikes, yspikes, 'Color', 'K');
            end 

            subplot(2,1,2); hold on
            title("Spike firing at soma")
            xlabel("Timestep")
 
            spikeTimes = find(nG.resultsAfter.didSomaSpikePerTimestep.' == 1);
            xspikes = repmat(spikeTimes, 3, 1);
            yspikes = nan(size(xspikes));
  
            if ~isempty(yspikes)
                   yspikes(1, :) = 0;
                   yspikes(2, :) = 1;
            end

            plot(xspikes, yspikes, 'Color', 'K');
        end
    end
end