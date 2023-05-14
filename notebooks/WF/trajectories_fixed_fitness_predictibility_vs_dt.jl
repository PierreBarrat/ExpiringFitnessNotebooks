### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
begin
	using Revise
	using Pkg; Pkg.activate("../../")
	using Chain
	using Distributions
	using FrequencyTrajectories
	using Plots
	using WrightFisher
	using StatsBase
end

# ╔═╡ 6cf21eab-ddb6-4c91-9e6a-d925d61702c2


# ╔═╡ 9e201d74-9e35-46a1-ad42-716630168318


# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
	N = 1_000_000
	L = 150
	μ = 0
	md"**Constants**"
end

# ╔═╡ c75346e5-181b-437d-9210-163cd6fdd995
s = 0.3

# ╔═╡ 2cc6dc91-c87a-4eb7-8707-0a7df2be6d1a
ρ = 1/4

# ╔═╡ b7f8bdbb-7c61-477f-a0eb-813b16701aac
md"## High fitness effect, high sweep rate"

# ╔═╡ 9009ba2e-1661-4ede-b01a-e29a71763e71
md"""
Why do trajectories become less predictible when the sampling time $\Delta t$ becomes larger? Let's choose one frequency bin and look at all rising trajectories entering this bin. 
I then plot two things: 
1. the average fitness of trajectories at the point where they enter the bin. 
2. for all trajectories entering the bin, the number that will go up at the next time point vs the number that will go down. 
"""

# ╔═╡ 86ff1f8d-b2e3-4a26-9e99-234cc4bfc38c
always_below = false

# ╔═╡ d82545fc-e789-4be5-a4f2-29354bf67244
md"""
The fitness of trajectories w.r. to the mean fitness of the population decreases as $\Delta t$ increases. This is expected as the very fit trajectories rise quickly, and are likely to not be seen in the bin for a high $\Delta t$. 
Inversely, we see many low fitness trajectories for a high value of $\Delta t$.
"""

# ╔═╡ e2daf6a9-d4de-4c65-8285-fd3896b9ed66
md"""
Below the number of trajectories going up vs going down after the first time they're seen in the bin. At higher $\Delta t$: 
- less go up, which is expected (fast growing ones are missed)
- more go down! This is very likely due to the ones turning around at a value higher than the bin, and caught on the way down for a high $\Delta t$. They are counted as going up when $\Delta t=1$. 
"""

# ╔═╡ e4943794-9de3-4564-967f-c71f29007965
md"# Functions"

# ╔═╡ a89608e8-511b-4217-b90d-7aa8de5ef69b


# ╔═╡ 8743fc5a-77d0-4975-8b9f-a607aa4387ec
begin
	function map_genotype_fitness!(T::FT.Trajectory, time_to_pop::Dict)
		T.ϕgenotype = map(T.t) do t
			N = sum(x -> x.count, time_to_pop[t])
			av_fit = sum(x -> x.fitness*x.count, time_to_pop[t])/N
	
			Nmut = 0
			av_fit_mut = 0
			for x in time_to_pop[t]
				if x.genotype.seq[T.pos] == T.state
					av_fit_mut += x.fitness*x.count
					Nmut += x.count
				end
			end
			av_fit_mut /= Nmut
			av_fit_mut - av_fit
			av_fit_mut - av_fit
		end
	end
	"""
		map_genotype_fitness!(trajectories, cb)

	For each trajectory: use the `gen_fitness` callback value to compute the average fitness of all genotypes carrying the corresponding mutation w.r. to the average fitness of the population. 
	"""
	function map_genotype_fitness!(trajectories, cb)
		time_to_pop = Dict(x.t => x.gen_fitness for x in cb)
		for T in trajectories
			map_genotype_fitness!(T, time_to_pop)
		end
		return nothing
	end
end

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(s, ρ; Δt = 1, T = max(500/ρ, 500/s))
	# setting parameters
	fitness_distribution = Exponential(s/2)	
	switchgen = round(Int, 1/ρ)

	# initial population
	H = rand(fitness_distribution, L)
	ϕ = AdditiveFitness(; L, H)
	pop = Pop(ϕ; N, L, μ)

	genotype_fitness(pop) = map(collect(pairs(pop.genotypes))) do (id, g)
		(
			genotype = g, 
			count = pop.counts[id],
			fitness = WF.fitness(g, pop),
		)
	end
	
	cb = (
		f1 = WF.frequencies,
		mut_fitness = WF.fields,
		gen_fitness = genotype_fitness,
		varpos_strict = pop -> WF.diversity(
			pop; method = :variable_positions, variable=1/N
		)/L,
		varpos_5 = pop -> WF.diversity(
			pop; method = :variable_positions, variable=0.05,
		)/L,
		n_genotypes = pop -> length(pop.genotypes),
	)
	
	cb_vals, switch_times = WF.Tools.evolve_sample!(
		pop, T, Δt, cb;
		fitness_distribution,
		switchgen, 
		change_init_field=true, 
		change_field_time = :periodic,
		max_freq=0., 
	)

	diversity = map(cb_vals) do x
		(varpos_strict = x.varpos_strict, varpos_5 = x.varpos_5)
	end

	f1 = zeros(Float64, length(cb_vals), 2*L)
	for (i, x) in enumerate(cb_vals)
		f1[i, :] .= x.f1
	end
	h = zeros(Float64, length(cb_vals), 2*L)
	for (i, x) in enumerate(cb_vals)
		h[i, :] .= x.mut_fitness
	end
	
	trajectories = get_trajectories(f1; time_scaling=Δt, fitness=h)
	map_genotype_fitness!(trajectories, cb_vals)
	
	sort!(trajectories, by = x -> x.t[1])
	sort!(switch_times, by = x -> x.t)

	return trajectories, diversity, switch_times
end

# ╔═╡ 832fe8a7-038d-48ca-94e1-64cf83c360c5
begin
	Δt = 4
	T, div, st = simulate(s, ρ; Δt, T = 5000)
	fb = FrequencyBin(0.5, 0.05)
	filter!(T, fb; always_below=true)

	
	p = plot()
	for t in T
		plot!(t, fb)
	end
	av_T, tvals, _ = mean(T, fb; always_below=true)
	plot!(tvals, av_T, line=(:black, 3), label="")
end

# ╔═╡ a5410237-9acf-4b51-b1a2-3f76087ef60a
function bin_stat(fb, Δt)
	T, div, st = simulate(s, ρ; Δt, T = 25_000)
	filter!(T, fb; always_below)

	fitness_values = map(T) do x
		i = findfirst(==(x.time_at_bin[fb]), x.t)
		x.ϕgenotype[i]
	end

	n_down = count(T) do x
		i = findfirst(==(x.time_at_bin[fb]), x.t)
		x.final_state != :missing && (x.f[i+1] < x.f[i])
	end

	n_up = count(T) do x
		i = findfirst(==(x.time_at_bin[fb]), x.t)
		x.final_state != :missing && (x.f[i+1] > x.f[i])
	end

	return n_up, n_down, fitness_values
end

# ╔═╡ 70025677-825f-46fb-8ba6-bdff6dc78a43
let
	fb = FrequencyBin(0.5, 0.05)
	p = plot(
		xlabel = "fitness", 
		title = "Fitness of trajectories at freq $(fb.f)", 
		frame=:box
	)
	for Δt in [1, 2, 5, 10, 20]
		_, _, fitness_values = bin_stat(fb, Δt)
		plot!(
			sort(fitness_values), (1:length(fitness_values))/length(fitness_values); label="Δt=$(Δt)"
		)
	end
	p
end

# ╔═╡ 3ad4ea81-e180-4ba4-a427-405c772cba6d
let
	fb = FrequencyBin(0.5, 0.05)
	p = plot(
		xlabel = "# down", 
		ylabel = "# up", 
		frame=:box
	)
	for Δt in [1, 2, 5, 10, 20]
		nup, ndown, fitness_values = bin_stat(fb, Δt)
		scatter!([ndown], [nup], label="$Δt")
	end
	p
end

# ╔═╡ 344605b7-7735-49e2-b965-37dac463ab85
md"# Tests"

# ╔═╡ ab7f8ac2-720b-4c76-8db2-1118ab753ebf


# ╔═╡ 848b216e-3759-4e79-a6f3-63163c58a3fe


# ╔═╡ ec2f318a-082c-4213-8597-a85e491071e3


# ╔═╡ 2fda3321-ff5a-4134-8646-7da8b2666c6e


# ╔═╡ Cell order:
# ╠═10b7789c-e56f-4b71-8e38-312af11bb751
# ╟─6cf21eab-ddb6-4c91-9e6a-d925d61702c2
# ╠═9e201d74-9e35-46a1-ad42-716630168318
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╠═c75346e5-181b-437d-9210-163cd6fdd995
# ╠═2cc6dc91-c87a-4eb7-8707-0a7df2be6d1a
# ╟─b7f8bdbb-7c61-477f-a0eb-813b16701aac
# ╠═832fe8a7-038d-48ca-94e1-64cf83c360c5
# ╟─9009ba2e-1661-4ede-b01a-e29a71763e71
# ╠═86ff1f8d-b2e3-4a26-9e99-234cc4bfc38c
# ╟─d82545fc-e789-4be5-a4f2-29354bf67244
# ╠═70025677-825f-46fb-8ba6-bdff6dc78a43
# ╟─e2daf6a9-d4de-4c65-8285-fd3896b9ed66
# ╠═3ad4ea81-e180-4ba4-a427-405c772cba6d
# ╠═e4943794-9de3-4564-967f-c71f29007965
# ╠═a5410237-9acf-4b51-b1a2-3f76087ef60a
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═a89608e8-511b-4217-b90d-7aa8de5ef69b
# ╠═8743fc5a-77d0-4975-8b9f-a607aa4387ec
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═ab7f8ac2-720b-4c76-8db2-1118ab753ebf
# ╠═848b216e-3759-4e79-a6f3-63163c58a3fe
# ╠═ec2f318a-082c-4213-8597-a85e491071e3
# ╠═2fda3321-ff5a-4134-8646-7da8b2666c6e
