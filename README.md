1/27 TODO
- change plasticity flag (1, 2, 3, 4)
- why is orientation tuning = 1?
- 


things I fixed:
- 
- 

need to do:
- account for > 4 orientations
- account for noisy L4 inputs
- documentation for methods
- write GET methods instead of making stim/dend instance properties public
- calling method on an instance? also is there an easier way to reference instance vars?
- write test cases for calcDendActivity / Dendritic_Clusters
- EASY: GROUP ALL RETURN VARS TOGETHER

QUESTIONS:
* NeuronGraphics --> layout2
- how was the heatmap colored? how were the blue bars made
- something wrong with histogram

* NeuronGraphics --> layout5
- what does 'branchWMax' mean? I only have wMax for the synapses

* ModelNeuron --> dendriticRFAnalyze*
- difference between Size1 and Size2?

* makeScan/do_Bar_Scan *
- adding noise: 
    - shouldn't this only be positive noise? you add 0.9 * enoise(1) (+/- 0.9)
    - do we want to add the same amount of noise to each L4? i fixed my code so it doesn't

* runModel/do_Run_BarStim *
- how were the hardcoded constants in stimulusUpdate (do_Run_BarStim) chosen?
- pls explain branchGNa calculation

* ModelNeuron --> CalcDendActivity *
- why do we add multiplicative noise to synapses? On what interval is this noise chosen? (enoise was -1,1) 
- potential bug: indexing record arrays at iTime, instead of recording timestep index

* ModelNeuron --> dendriticRFAnalyze
- what is a tuning curve?

PROGRESS:
- wrote and debugged makeStimulus()
- 

FUTURE CHANGES:
- examine effects of noise
- effects of unequal number of synapses per dendrite