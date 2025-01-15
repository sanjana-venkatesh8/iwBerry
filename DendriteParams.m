classdef DendriteParams
    %DENDRITEPARAMS Stores dendrite parameters
    % Dendrite parameters:
    %     nSynapses - (int) # of synapses
    %     nDendrites - (int) # of basal dendrites
    %     branchSize - (int) # of synapses per branch
    %     synStrength - (?) max synapse strength (mV)
    %     synAttenuation - (?) attenuation at soma
    %     branchThresh - (?) threshold for NMDA spike (mV)
    %     branchStrength
    %     somaThresh
    %     backpropAP
    %     synNoise
    %     branchNoise
    %     somaNoise
    %     weightsRange
    %     duPotent
    %     duDepress
    %     duDecay
    %     duBaseline
    %     uRecycle
    %     uMax
    %     scaleNMDA
    %     scaleNoNMDA
    %     branchGLeak % TODO: what are the g_s? conductances
    %     initSynGInhib - initial inhibitory conductance of a synapse
    %     somaPlus
    %     somaMinus
    %     somaGLeak
    %     initSomaGInhib
    % TODO: fill in remaining property descriptions
    
    properties
        % input parameters
        nSynapses;
        nDendrites; % TODO: remove all occurrences of this
        branchSize;
        synStrength;
        synAttenuation;
        branchThresh;
        branchStrength;
        somaThresh;
        backpropAP;
        synNoise;
        branchNoise;
        somaNoise;
        weightsRange;
        duPotent;
        duDepress;
        duDecay;
        duBaseline;
        uRecycle;
        uMax;
        plastTime;
        scaleNMDA double {mustBeNonpositive} % TODO: is this a requirement?
        scaleNoNMDA double {mustBeNonnegative} % ^
        branchGLeak;
        initSynGInhib;
        somaPlus {mustBeNonnegative}
        somaMinus {mustBeNonpositive}
        somaGLeak;
        initSomaGInhib;

        % non-input parameters
        nBranches;
    end
    
    % TODO: remove 'non-leaky branch' conditionals
    methods (Access = public)
        function dendParams = DendriteParams(nSynapses, nDendrites, ...
            branchSize, synStrength, synAttenuation, branchThresh, ...
            branchStrength, somaThresh, backpropAP, synNoise, ...
            branchNoise, somaNoise, weightsRange, duPotent, duDepress, ...
            duDecay, duBaseline, uRecycle, uMax, plastTime, scaleNMDA, ...
            scaleNoNMDA, branchGLeak, initSynGInhib, somaPlus, ...
            somaMinus, somaGLeak, initSomaGInhib)
            
            %DENDRITEPARAMS Construct an instance of this class
            dendParams.nSynapses = nSynapses;
            dendParams.nDendrites = nDendrites;
            dendParams.branchSize = branchSize;
            dendParams.synStrength = synStrength;
            dendParams.synAttenuation = synAttenuation;
            dendParams.branchThresh = branchThresh;
            dendParams.branchStrength = branchStrength;
            dendParams.somaThresh = somaThresh;
            dendParams.backpropAP = backpropAP;
            dendParams.synNoise = synNoise;
            dendParams.branchNoise = branchNoise;
            dendParams.somaNoise = somaNoise;
            dendParams.weightsRange = weightsRange;
            dendParams.duPotent = duPotent;
            dendParams.duDepress = duDepress;
            dendParams.duDecay = duDecay;
            dendParams.duBaseline = duBaseline;
            dendParams.uRecycle = uRecycle;
            dendParams.uMax = uMax;
            dendParams.plastTime = plastTime;
            dendParams.scaleNMDA = scaleNMDA;
            dendParams.scaleNoNMDA = scaleNoNMDA;
            dendParams.branchGLeak = branchGLeak;
            dendParams.initSynGInhib = initSynGInhib;
            dendParams.somaPlus = somaPlus;
            dendParams.somaMinus = somaMinus;
            dendParams.somaGLeak = somaGLeak;
            dendParams.initSomaGInhib = initSomaGInhib;
            
            dendParams.nBranches = ceil(nSynapses / branchSize);
            % TODO: what to do if missing params in input?? make name-value
            % args
        end

        % TODO: insert functions here if needed
    end
end

