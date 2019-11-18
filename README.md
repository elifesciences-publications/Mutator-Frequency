This repository contains the associated simulation code for the manuscript “Selection on mutators is not frequency-dependent” 
by Yevgeniy Raynes and Daniel M. Weinreich.https://doi.org/10.7554/eLife.51177

Simulation code is based on the individual-based, stochastic simulation previously described in Raynes et al. (PNAS 115(13): 3422-3427). 
In brief, we consider haploid asexual populations of constant size, N, evolving in discrete, non-overlapping generations according to 
the Wright-Fisher model (Ewens, 2004). Populations are composed of genetic lineages - subpopulations of individuals with the same genotype 
(stored in mutable struct Lineage). A genotype is modelled as an array of 99 fitness-affecting loci and 1 mutation rate modifier locus, 
which in a mutator state raises the genomic mutation rate m-fold (defined by mutator_strength). We assume constant fitness effects: 
beneficial mutations at the fitness loci increase fitness by a constant effect sb, while deleterious mutations decrease fitness by a 
constant effect sd. We assume additive fitness effects and so calculate fitness of a lineage as the sum of fitness effects of all of 
its mutations). Simulations start with the mutator allele at a frequency of x_0 (defined by init_mut_N) and continue until it either 
fixes (reaches the frequency of 1.0) or is lost (reaches the frequency of 0.0) from a population.

Every generation the size of each lineage i is randomly sampled from a multinomial distribution with expectation given by the product of 
its size and relative fitness. Reproduction is encoded by the functon wright_fisher_reproduction. Upon reproduction, each lineage acquires 
a random number of fitness-affecting mutations M, drawn from a Poisson distribution with mean equal to the product of its size and its 
total per-individual mutation rate, (Ub+Ud), where Ub and Ud are the deleterious and beneficial mutation rates respectively. 
The number of beneficial and deleterious mutationsis then drawn from a binomial distribution with n=M and P = Ub⁄((Ub+Ud)) and 
new mutations are assigned to randomly chosen non-mutated fitness loci. Mutation is encoded by the function mutate_population.

Parameter values for a given simulation run are defined in Lines 117 – 123 as follows: 

	pop_Ni = 10e7 			(Total population size)
	init_mut_N = n0 		(Size of the initial mutator subpopulation) 
	sb = 0.1 			(Selective effect of new beneficial mutations, constant)
	sd = 0.1 			(Selective effect of new deleterious mutations, constant)
	mutator_strength = 100.0 	(The fold increase in the mutator mutation rate over the non-mutator)
	Ub = 0.000001 			(Per-individual beneficial mutation rate)
	Ud = 0.0001 			(Per-individual deleterious mutation rate)

Simulation is run in Julia 1.0 and produces two types of files (where X is the number of the particular simulation run in a job array if needed):
	
	Xmutator_traces.csv: records mutator frequency for every generation of the simulation 
                       (comma-delimited, different runs of the simulation are separated by line breaks).
	
	Xtime_to_mut.csv: records the time to mutator fixation (0.0 if mutator is lost) in each simulation run 
                    (comma-delimited, sets of simulation runs at different starting frequencies are separated by line breaks).
