classdef ModelResults
    %ModelResults Stores results after running Jens2Plasticity
    
    properties
        nTotalSpikes
        NMDASpikesPerTimestep
        somaActivityPerTimestep
        somaEPSPsPerTimestep
        didSomaSpikePerTimestep
        branchNMDASpikeRate
        avgVoltageRespPerBranch
        branchHistogram
        somaHistogram
        branchRF
        cumulStimRF
        somaThreshPerTimestep
        somaGainPerTimestep
        somaInhibPerTimestep
        synURecord
        synL4Record
        branchWMaxRecord
        branchVRecord
        branchSpikeRecord
    end
    
    methods
        function resultsObj = ModelResults(nTotalSpikes, NMDASpikesPerTimestep, ...
                somaActivityPerTimestep, somaEPSPsPerTimestep, didSomaSpikePerTimestep, ...
                branchNMDASpikeRate, avgVoltageRespPerBranch, ...
                branchHistogram, somaHistogram, branchRF, cumulStimRF, ...
                somaThreshPerTimestep, somaGainPerTimestep, somaInhibPerTimestep, ...
                synURecord, synL4Record, branchWMaxRecord, branchVRecord, branchSpikeRecord)
            %ModelNeuron Construct an instance of this class
            resultsObj.nTotalSpikes = nTotalSpikes;
            resultsObj.NMDASpikesPerTimestep = NMDASpikesPerTimestep;
            resultsObj.somaActivityPerTimestep = somaActivityPerTimestep;
            resultsObj.somaEPSPsPerTimestep = somaEPSPsPerTimestep;
            resultsObj.didSomaSpikePerTimestep = didSomaSpikePerTimestep;
            resultsObj.branchNMDASpikeRate = branchNMDASpikeRate;
            resultsObj.avgVoltageRespPerBranch = avgVoltageRespPerBranch;
            resultsObj.branchHistogram = branchHistogram;
            resultsObj.somaHistogram = somaHistogram;
            resultsObj.branchRF = branchRF;
            resultsObj.cumulStimRF = cumulStimRF;
            resultsObj.somaThreshPerTimestep = somaThreshPerTimestep;
            resultsObj.somaGainPerTimestep = somaGainPerTimestep;
            resultsObj.somaInhibPerTimestep = somaInhibPerTimestep;
            resultsObj.synURecord = synURecord;
            resultsObj.synL4Record = synL4Record;
            resultsObj.branchWMaxRecord = branchWMaxRecord;
            resultsObj.branchVRecord = branchVRecord;
            resultsObj.branchSpikeRecord = branchSpikeRecord;
        end
    end
end

