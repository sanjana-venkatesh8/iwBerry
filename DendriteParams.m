classdef DendriteParams
    %DENDRITEPARAMS Stores dendrite parameters
    
    properties
        % input parameters
        nSynapses;                                              % number of synapses
        branchSize;                                             % number of synapses per branch
        synStrength;                                            % max synapse strength (mV) - TODO: REMOVE; this isn't being used
        synAttenuation;                                         % attenuation at soma (mV)
        branchThresh;                                           % threshold for an NMDA spike (mV)
        branchStrength;                                         % strength of branch NMDA spike at soma (mV)
        somaThresh;                                             % threshold for a somatic spike (mV)
        backpropAP;                                             % amplitude of backpropagating action potential (mV)
        synNoise double {mustBeInRange(synNoise, 0, 1)}         % noise on synapses
        branchNoise double {mustBeInRange(branchNoise, 0, 1)}   % noise on branches
        somaNoise double {mustBeInRange(somaNoise, 0, 1)}       % noise at soma
        weightsRange;                                           % range for initial weights
        duPotent;                                               % change in internal potential (u) for   YES presynaptic and YES postsynaptic spike
        duDepress;                                              % change in internal potential (u) for  YES presynaptic and NO  postsynaptic spike
        duDecay;                                                % change in internal potential (u) for    NO presynaptic  and YES postsynaptic spike
        duBaseline;                                             % change in internal potential (u) for NO presynaptic  and NO  postsynaptic spike
        uRecycle;                                               % threshold internal potential (u) for a spine recycle
        uMax;                                                   % maximum internal potential (u)
        plastTime;                                              % size of plasticity time window when potentiated
        scaleNMDA % double {mustBeNonpositive}                    % synaptic scaling when there is an NMDA spike
        scaleNoNMDA % double {mustBeNonnegative}                  % synaptic scaling when no NMDA spike
        branchGLeak;                                            % direct modulation of weights (0) or inhibitory shunting-based control (1)
        initSynGInhib;                                          % initial inhibitory conductance (G) of synapse
        somaPlus {mustBeNonnegative}                            % somatic threshold change if there is a somatic spike
        somaMinus {mustBeNonpositive}                           % somatic threshold change if no somatic spike
        somaGLeak;                                              % leak conductance at soma
        initSomaGInhib;                                         % initial inhibitory conductance (G) of soma)

        % non-input parameters
        nBranches;                                              % number of branches in neuron
    end
    
    % TODO: remove 'non-leaky branch' conditionals
    methods (Access = public)
        function dendParams = DendriteParams(args)
            %DENDRITEPARAMS Construct an instance of this class

            arguments
                args.nSynapses;
                args.branchSize;
                args.synStrength;
                args.synAttenuation;
                args.branchThresh;
                args.branchStrength;
                args.somaThresh;
                args.backpropAP;
                args.synNoise;
                args.branchNoise;
                args.somaNoise;
                args.weightsRange;
                args.duPotent;
                args.duDepress;
                args.duDecay;
                args.duBaseline;
                args.uRecycle;
                args.uMax;
                args.plastTime;
                args.scaleNMDA
                args.scaleNoNMDA
                args.branchGLeak;
                args.initSynGInhib;
                args.somaPlus
                args.somaMinus
                args.somaGLeak;
                args.initSomaGInhib;
            end

            dendParams.nSynapses = args.nSynapses;
            dendParams.branchSize = args.branchSize;
            dendParams.synStrength = args.synStrength;
            dendParams.synAttenuation = args.synAttenuation;
            dendParams.branchThresh = args.branchThresh;
            dendParams.branchStrength = args.branchStrength;
            dendParams.somaThresh = args.somaThresh;
            dendParams.backpropAP = args.backpropAP;
            dendParams.synNoise = args.synNoise;
            dendParams.branchNoise = args.branchNoise;
            dendParams.somaNoise = args.somaNoise;
            dendParams.weightsRange = args.weightsRange;
            dendParams.duPotent = args.duPotent;
            dendParams.duDepress = args.duDepress;
            dendParams.duDecay = args.duDecay;
            dendParams.duBaseline = args.duBaseline;
            dendParams.uRecycle = args.uRecycle;
            dendParams.uMax = args.uMax;
            dendParams.plastTime = args.plastTime;
            dendParams.scaleNMDA = args.scaleNMDA;
            dendParams.scaleNoNMDA = args.scaleNoNMDA;
            dendParams.branchGLeak = args.branchGLeak;
            dendParams.initSynGInhib = args.initSynGInhib;
            dendParams.somaPlus = args.somaPlus;
            dendParams.somaMinus = args.somaMinus;
            dendParams.somaGLeak = args.somaGLeak;
            dendParams.initSomaGInhib = args.initSomaGInhib;
            
            dendParams.nBranches = ceil(args.nSynapses / args.branchSize);
        end
    end
end

