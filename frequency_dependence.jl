using Distributions

mutable struct Lineage
  fitness::Float64
  size::Int64
  state::Vector{Float64}
  mut_rate::Float64
end

function average_fitness(population::Dict{Array{Float64,1}, Lineage})
#calculate average fitness of the population. also outputs population size
  popN = Int64
  popN = 0
  popW = Float64
  popW = 0
  for line in values(population)
    popN+=line.size
    popW+=line.size * line.fitness
  end
  popW = popW/popN
  return popW, popN
end

function assay_mutation_rate(population::Dict{Array{Float64,1}, Lineage})
#calculates mutator frequency
  popN = Int64
  popN = 0
  mutN = Int64
  mutN = 0
  for line in values(population)
    popN+=line.size
    if line.mut_rate>1 mutN+=line.size end
  end
  mutF = Float64
  mutF = mutN/popN
  return mutF
end

function avail_sites(state::Array{Float64,1})
  sites = Int64[]
  for i = 1:length(state)
    if state[i] == 0 push!(sites,i) end
  end
  return sites
end

function mutate_population(population::Dict{Array{Float64,1}, Lineage}, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64, mut_multiplier::Float64)
  new_population = Dict{Array{Float64,1}, Lineage}()
  #initialize new population dictionary
  for (key, line) in population
    mutations = rand(Poisson(line.mut_rate * (Ud + Ub) * line.size))
    if mutations > line.size mutations = line.size end
    bmutations = rand(Binomial(mutations, Ub/(Ub+Ud)))
    dmutations = mutations - bmutations

    open_sites =   avail_sites(line.state)
    #generates a list of indexes of all loci that are available for mutation
    mutation_positions = open_sites[rand(1:end,(dmutations+bmutations))] #generates positions of new mutations by randomly picking indexes from open_sites

      #for each new mutation, copies the state, adds mutation and generates the key.
      #check if it is alrady in the new_population. if it is, add 1 individual to it, if it's not - make another lineage with
      #the key and size 1
    c=1
    for i = 1:bmutations
      #assign beneficial mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] = sb
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end

    for i = 1:dmutations
      #assign deleterious mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] = -sd
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end

    new_size = line.size - bmutations - dmutations
    if  new_size > 0
      if key in keys(new_population) new_population[key].size+=new_size
      else new_population[key] = Lineage(line.fitness, new_size, line.state, line.mut_rate)
      end
    end
  end
  return new_population
end

function wright_fisher_reproduction(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[] #initialize lineage liste
  for (key,line) in population
    push!(proby_list, max(line.size/N0 * line.fitness/popw, 0.0))
    #representation of a lineage in the next generations depends on its frequency and relative fitness
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mut_rate)
    end
  end
  return new_population
end

function simulate(n0,job_id)
  pop_Ni = 10000000
  init_mut_N = n0
  sb = 0.1
  sd = 0.1
  mutator_strength = 100.0
  Ub = 0.000001
  Ud = 0.0001

  population = Dict{Array{Float64,1}, Lineage}()
  wt_state = zeros(Float64,100)
  wt_state[1] = 1.0
  population[wt_state] = Lineage(1,pop_Ni-init_mut_N,wt_state,1.0)
  mut_state = copy(wt_state)
  mut_state[1] = mutator_strength
  population[mut_state] = Lineage(1,init_mut_N,mut_state,mutator_strength)
  mut_f = assay_mutation_rate(population)
  popw, popN = average_fitness(population)
  generations = Int64
  generations = 0
  #initialize generations counter
  mut_trajectory = Float64[]
  push!(mut_trajectory, mut_f)

  while 0.0<mut_f<1.0
    popw, pop_N = average_fitness(population)
    population = wright_fisher_reproduction(population, pop_N, pop_Ni, popw)
    population = mutate_population(population, Ub, sb, Ud, sd, mutator_strength)
    mut_f = assay_mutation_rate(population)
    popw, pop_N = average_fitness(population)
    push!(mut_trajectory, mut_f)
    generations+=1
    #println("At the end of generation: ", generations, " Mut frequency: ", mut_f, " Fitness: ", popw, " Population size: ", pop_N, " Number of lines ", length(population))

  end
  outfile = open(string(job_id,"mutator_traces.csv"), "a")
  write(outfile, join(mut_trajectory, ","), "\n")
  close(outfile)
  return mut_f ==1.0, generations, popw
end

job_id = ARGS[1]
for n0 in [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
  time_to_mut= Float64[]
  #time_to_mut list records all times of successful mutator hitchhiking
  for run = 1:25
    output = simulate(n0, job_id)
    if output[1]
      push!(time_to_mut, output[2])
    else
      push!(time_to_mut, 0)
    end
  end
  outfile = open(string(job_id,"time_to_mut.csv"), "a")
  write(outfile, join(time_to_mut, ","), "\n")
  close(outfile)
end
