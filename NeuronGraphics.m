classdef NeuronGraphics
    properties
        resultsBefore ModelResults
        resultsPlast ModelResults
        resultsAfter ModelResults
        
        resultsRFBefore RFResults
        resultsRFAfter RFResults
    end

    methods
        function neuronGraphics = NeuronGraphics(...
                resultsBefore, resultsPlast, resultsAfter, resultsRFBefore, resultsRFAfter)
            arguments
                resultsBefore ModelResults
                resultsPlast ModelResults
                resultsAfter ModelResults
                
                resultsRFBefore RFResults
                resultsRFAfter RFResults
            end

            neuronGraphics.resultsBefore = resultsBefore;
            neuronGraphics.resultsPlast = resultsPlast;
            neuronGraphics.resultsAfter = resultsAfter;
            
            neuronGraphics.resultsRFBefore = resultsRFBefore;
            neuronGraphics.resultsRFAfter = resultsRFAfter;
        end

        function layout1(nG)
        %LAYOUT1 Plots orientation tuning vs. NMDA spike rate, colored by 
        % RF size, before and after plasticity.
            % TODO: might need * 100 for branchIOrient
            figure(Name="Layout 1: Orientation Tuning Across Subunits")            
            subplot(2, 1, 1)
            scatter(nG.resultsBefore.branchNMDASpikeRate, ...
                nG.resultsRFBefore.branchIOrient, 50, nG.resultsRFBefore.branchSize1, 'filled');
            c = colorbar;
            c.Label.String = "RF Size (pixels)";
            title("Before plasticity")
            ylabel("Orientation Tuning Index")
            xlabel("NMDA Spike Rate")

            subplot(2, 1, 2)
            scatter(nG.resultsAfter.branchNMDASpikeRate, ...
                nG.resultsRFAfter.branchIOrient, 50, nG.resultsRFAfter.branchSize1, 'filled');
            c = colorbar;
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
                                "Orientation index = %d\n" + ...
                                "Receptive field size = %1.1d\n" + ...
                                "Max EPSP = %d mV\n" + ...
                                "Max response = %d", ...
                                iBranch, nG.resultsRFAfter.branchIOrient(iBranch), ...
                                nG.resultsRFAfter.branchSize1(iBranch), NaN, NaN);
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

            % TODO: check if plotting right hist.
            % plot histogram of ???
            nexttile([1 2])
            nOrient = size(nG.resultsRFBefore.allBranchOrient, 2);
            histogram(nG.resultsRFAfter.allBranchOrient(iBranch, :), BinEdges=0.5+(0:1:nOrient))
            xticks(0:1:nOrient)
            title("Branch orientation tuning")
            xlabel("Bar orientation")

            % use compositeOrient and compositeTuning for histogram
            % heat map
        end

        function layout3(nG)
        %LAYOUT3 Plots the soma receptive field vs. composite receptive
        %field
            figure(Name="Layout 3: Soma Receptive Field vs. Composite Receptive Field")
            
        end
        
        function layout4(nG, iBranch, branchGLeak, tau)
        %LAYOUT4 Plots neuron adaptation: NMDA rate, synaptic max weight, and branch
        %inhibition over time

        % tau = smoothing window
        % TODO: rateTarget = NMDA firing rate target

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
        figure(Name=title4)
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

            % Calculate synaptic weights
            branchExc = 1 ./ (1 + exp(-1 .* nG.resultsPlast.synURecord(:, iSyn, iBranch)));
            % TODO: branchInh might be wrong
            % TODO: where is 60 from?
            branchInh = 60 ./ nG.resultsPlast.branchWMaxRecord(iBranch) - 1 - branchGLeak;
            synWeights = 60 * branchExc / (branchExc + branchInh + branchGLeak);

            title5 = sprintf("Layout 5: Plasticity in Synapse %d (x = %d, y = %d)", iSyn, NaN, NaN);
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
            % histogram(nG.resultsBefore.)

            % UHistogram_Before = histcounts(modelNeuron.synUInput, NumBins=51);
            % WMaxHistogram_Before = histcounts(modelNeuron.synInputWMax, NumBins=21);
        end
    end
end