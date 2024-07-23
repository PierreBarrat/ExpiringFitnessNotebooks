### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 58d01308-3dc9-11ef-0b3a-99e848c2ffd2
begin
	using Revise
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using WrightFisher
	using StatsBase

	using Logging, LoggingExtras
end

# ╔═╡ b0aaad98-3c42-4693-971c-6c4f4512364e


# ╔═╡ f1e3f20b-9975-4434-80f4-883fab8e74f6
notebook_name = replace(@__FILE__, r"\.jl#==#.*" => "") |> basename

# ╔═╡ 6bf1d768-e2d0-4f15-9d51-8582f4a27879
savedir = "data_" * notebook_name * "/"

# ╔═╡ c02edeec-3547-43d5-b132-a8039205488f
mkpath(savedir)

# ╔═╡ 5fb29a1a-3886-41e5-96c4-59d6c10d3594
begin
	N = 100_000
	L = 200
	μ = 0
	f0 = 0.01
	md"**Constants**"
end

# ╔═╡ bc0a2f1a-e9e1-4716-947a-bd1d88a5a524
function simulate(s, ρ; Δt = 1, change_field_time = :periodic)
	# setting parameters
	fitness_distribution = Exponential(s/2)

	T = max(1000/ρ, 1000/s)
	
	switchgen = round(Int, 1/ρ)

	# initial population
	H = rand(fitness_distribution, L)
	ϕ = AdditiveFitness(; L, H)
	pop = Pop(ϕ; N, L, μ)
	cb = (
		f1 = WF.frequencies,
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
		change_field_time,
		f0 = f0,
		max_freq = f0,
	)

	diversity = map(cb_vals) do x
		(varpos_strict = x.varpos_strict, varpos_5 = x.varpos_5)
	end

	f1 = zeros(Float64, length(cb_vals), 2*L)
	for (i, x) in enumerate(cb_vals)
		f1[i, :] .= x.f1
	end
	trajectories = get_trajectories(f1; time_scaling=Δt) |> FrequencyTrajectories.trajectories_to_dataframe

	return trajectories, diversity, switch_times
end

# ╔═╡ e61a5c4e-efa4-4c98-9a12-dcd9bbd283d0
md"## Utilities"

# ╔═╡ a1b0cb8f-2dce-473a-b6d3-70bb46612685


# ╔═╡ 475dfdd1-c9bd-4c5a-b6b7-159f7aa7c734
logrange(x, y, L) = map(exp, range(log(x), log(y), length=L))

# ╔═╡ 0d0ae104-9aaf-4419-a604-1996a4ab0c78
begin
	svals = [.02]
	ρvals = logrange(1e-1, 10*sqrt(10), 5) * svals[1]
	# ρvals = [10.] * svals[1]
	Δtvals = [10]
	cf_style = [:periodic, :random]
	md"**Variable parameters**"
end

# ╔═╡ dce3af86-e2bb-42a1-adb9-d56ae2d1aa16
parameters=map(Iterators.product(Δtvals, ρvals, svals, cf_style)) do (Δt, ρ, s, cfs)
	(s=s, ρ=ρ, Δt=Δt, cfs = cfs)
end |> (x -> vcat(x...));

# ╔═╡ 05acebd0-6f54-44a6-a308-37b73feee7df
trajectories, diversity = let
	tjs = Dict()
	div = Dict()
	for (i, p) in enumerate(parameters)
		@info p i/length(parameters)
		@time T, d = simulate(p.s, p.ρ; Δt=p.Δt, change_field_time = p.cfs)
		tjs[p] = T
		
		df = DataFrame()
		for x in d
			push!(df, x)
		end
		div[p] = df
	end
	tjs, div
end;

# ╔═╡ e317fe0b-a9d0-4cbc-9b18-cf5c9e816ba5
collect(keys(trajectories))

# ╔═╡ 5983190e-250e-4f3a-a563-0e6b73228202
filenames = map(enumerate(collect(keys(trajectories)))) do (idx, p)
	idx => (
		ρ = p.ρ,
		s = p.s,
		Δt = p.Δt,
		cfs = p.cfs,
		L=L, 
		N=N,
		trajectory_file = "trajectory_$(idx).csv",
		diversity_file = "diversity_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ ae7f50d2-0c73-4448-b9fd-289b27d9b07f
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, Dict(filenames))
end

# ╔═╡ 8e393026-a133-4314-988d-50c7d8f18664
for (i, v) in filenames
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[v.key])
	end
	open(savedir * v.trajectory_file, "w") do io
		CSV.write(io, trajectories[v.key])
	end
end

# ╔═╡ 3c96b374-a697-42e2-bc21-aaecd265d9ac
logrange(1e-1, 100, 4)

# ╔═╡ Cell order:
# ╠═58d01308-3dc9-11ef-0b3a-99e848c2ffd2
# ╠═b0aaad98-3c42-4693-971c-6c4f4512364e
# ╠═f1e3f20b-9975-4434-80f4-883fab8e74f6
# ╠═c02edeec-3547-43d5-b132-a8039205488f
# ╠═6bf1d768-e2d0-4f15-9d51-8582f4a27879
# ╠═5fb29a1a-3886-41e5-96c4-59d6c10d3594
# ╠═0d0ae104-9aaf-4419-a604-1996a4ab0c78
# ╠═3c96b374-a697-42e2-bc21-aaecd265d9ac
# ╠═dce3af86-e2bb-42a1-adb9-d56ae2d1aa16
# ╠═bc0a2f1a-e9e1-4716-947a-bd1d88a5a524
# ╠═05acebd0-6f54-44a6-a308-37b73feee7df
# ╠═e317fe0b-a9d0-4cbc-9b18-cf5c9e816ba5
# ╠═5983190e-250e-4f3a-a563-0e6b73228202
# ╠═ae7f50d2-0c73-4448-b9fd-289b27d9b07f
# ╠═8e393026-a133-4314-988d-50c7d8f18664
# ╟─e61a5c4e-efa4-4c98-9a12-dcd9bbd283d0
# ╠═a1b0cb8f-2dce-473a-b6d3-70bb46612685
# ╠═475dfdd1-c9bd-4c5a-b6b7-159f7aa7c734
