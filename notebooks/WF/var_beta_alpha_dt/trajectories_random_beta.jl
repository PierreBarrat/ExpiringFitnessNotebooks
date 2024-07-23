### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
begin
	using Revise
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using DelimitedFiles
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using WrightFisher
	using StatsBase
end

# ╔═╡ 19759f0e-ced3-11ed-2b40-8b01992454e2
notebook_name = replace(@__FILE__, r"\.jl#==#.*" => "") |> basename

# ╔═╡ d8666482-d8a4-43f4-ad21-81f6e972ab15
savedir = "data_" * notebook_name * "/"

# ╔═╡ 4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
mkpath(savedir)

# ╔═╡ 6cf21eab-ddb6-4c91-9e6a-d925d61702c2
md"""
## Idea

We keep constant $N_e$, $\langle \beta^2 \rangle$ and $\rho = 1/N_e\langle \beta^2 \rangle$. The distribution of $\beta$ is then parametrized by its mean $\langle\beta\rangle$ (details in the `get_β_distribution` function). We vary three parameters
-  $\langle \beta^2 \rangle < \langle\beta\rangle < \sqrt{\langle \beta^2 \rangle}$ (`mβ` in the code)
-  $\alpha$
-  $\Delta t$
"""

# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
	N = 1_000_000
	L = 150
	μ = 0
	Ne = 100
	β2 = .1
	ρ = 1/Ne/β2
	md"**Constants**"
end

# ╔═╡ c7909693-203b-49a8-8713-387dc287e4a2
md"""
The duration of a trajectory should be of order $T \simeq 1/(\rho\langle\beta^2\rangle)$. The number of partial sweeps happening during this time is $1/\langle\beta^2\rangle$. For this reason we expect to always have a supply of mutations if $L \gg 1/\langle\beta^2\rangle$. 
"""

# ╔═╡ a37e9a1e-e06f-4402-abf3-f82e22b63fe3
md"""
I will try the following: 
- vary $\beta$ (constant) and $\rho$ to some set of values;
- fix other parameters;

and try to see whether the supply of mutation runs out in each case. 

For a given $\beta$ and $\rho$, the simulation time should be $\gg 1/\rho/\beta^2$.
"""

# ╔═╡ afe7bf04-d998-46fc-afed-8e2cd21361d2
begin
	mβvals = collect(range(β2, sqrt(β2), length=7))[2:end-1]
	αvals = [.03, .1, .3]
	Δtvals = [1, 3, 10, 30]
	md"**Variable parameters**"
end

# ╔═╡ 1110d7fa-9c96-4b40-ae60-2396a08d4db0
parameters = map(Iterators.product(Δtvals, αvals, mβvals)) do (Δt, α, mβ)
	(mβ=mβ, α=α, Δt=Δt)
end |> (x -> vcat(x...));

# ╔═╡ f77e8aa9-7843-4917-a745-0fa3e6611939
"""
	get_β_distribution(mβ, β2)

Return a Beta distribution with mean `mβ` and second moment `β2`.
"""
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return Beta(a,b)
end

# ╔═╡ b7d9a8d0-476d-4d87-acd8-a06c7766198b
begin
	struct FitnessDistribution <: Sampleable{Univariate, Continuous}
		β :: Distribution
		α :: Float64
	end

	function FitnessDistribution(mβ, β2, α)
		@assert β2 < mβ < sqrt(β2) "Incorrect value <β>=$(mβ). The range of possible values is $(β2) < mβ < $(sqrt(β2))"

		return FitnessDistribution(get_β_distribution(mβ, β2), α)
	end
end

# ╔═╡ 880c9804-d0e3-4f71-a34c-fcef988e3953
begin
	import Base.rand
	function rand(rng::Distributions.AbstractRNG, x::FitnessDistribution)
		β = @chain rand(x.β) min(_, .999) max(_, .001) # for the log to be ok
		return -x.α * log(1-β)
	end
end

# ╔═╡ 420fceea-4103-4262-95f8-0fb4450da83f
let
	α = .1
	p = FitnessDistribution(.101, β2, α)
	S = rand(p, 1000)/α
end

# ╔═╡ 51b03b84-3946-4079-a43a-7ae0050f2b93
md"""
An important ratio is $s/\rho$: 
- if $s \gg \rho$, the sweeps happen one after the other;
- if $s/\rho \ll 1$, they overlap a lot. 

As a function of $\beta$, we have

$$s/\rho = -α\rho^{-1}\cdot\log(1-\beta).$$
"""

# ╔═╡ cb1a8bf0-776d-41d1-8a50-90c81292c5be
md"""

It could be interesting to see what happens in the regime $s/ρ < 1$ and $Lβ^2 \gg 1$, meaning that we have a large supply of mutations according but the sweeps overlap. 

In terms of $β$ and $ρ$ this regime is defined by 

$$\beta < 1 - e^{-\rho/\alpha} \;\;\text{and}\;\; \beta \gg L^{-1/2}$$

This is of course only possible if 

$$1 - e^{-\rho/\alpha} \gg L^{-1/2}$$
"""

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(mβ, α; Δt = 1)
	# setting parameters
	fitness_distribution = FitnessDistribution(mβ, β2, α)

	T = max(500/ρ, 100 * Ne)
	
	switchgen = round(Int, 1/ρ)

	# initial population
	H = zeros(Float64, L)
	ϕ = ExpiringFitness(; L, H, α)
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
		change_field_time = :random,
		max_freq=0., 
	)

	diversity = map(cb_vals) do x
		(varpos_strict = x.varpos_strict, varpos_5 = x.varpos_5)
	end

	f1 = zeros(Float64, length(cb_vals), 2*L)
	for (i, x) in enumerate(cb_vals)
		f1[i, :] .= x.f1
	end
	trajectories = get_trajectories(f1; time_scaling=Δt) |> FrequencyTrajectories.trajectories_to_dataframe

	return trajectories, diversity, f1, switch_times
end

# ╔═╡ e8a8acd2-0aa0-4c36-b7be-8ef356fc6692
trajectories, diversity, allele_frequencies, switch_info = let
	tjs = Dict()
	div = Dict()
	freq = Dict()
	switches = Dict()
	for (i, p) in enumerate(parameters)
		@info p i/length(parameters)
		@time T, d, f, s = simulate(p.mβ, p.α; Δt=p.Δt)
		# trajectories
		tjs[p] = T
		# diversity
		df = DataFrame()
		for x in d
			push!(df, x)
		end
		div[p] = df
		# frequencies
		freq[p] = f
		# switch info
		df = DataFrame()
		for x in s
			push!(df, x; promote=true)
		end
		for col in names(df)
			transform!(df, col => (X -> map(x -> isnothing(x) ? missing : x, X)) => col)
		end
		switches[p] = df
	end
	tjs, div, freq, switches
end;

# ╔═╡ 50075364-ebd7-43b8-b1b0-3aa0375b01b4
filenames = map(enumerate(collect(keys(trajectories)))) do (idx, p)
	idx => (
		Ne = Ne,
		β2 = β2,
		ρ = ρ,
		mβ = p.mβ,
		α = p.α,
		Δt = p.Δt,
		L=L, 
		N=N,
		trajectory_file = "trajectory_$(idx).csv",
		diversity_file = "diversity_$(idx).csv",
		allele_freq_file = "allele_frequencies_$(idx).dat",
		switch_info_file = "switch_info_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ 8172b776-f207-41fe-bd85-8fd2fada723f
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, Dict(filenames))
end

# ╔═╡ faf656e1-a817-4a30-a9e8-f86cb37ef35d
for (i, v) in filenames
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[v.key])
	end
	open(savedir * v.trajectory_file, "w") do io
		CSV.write(io, trajectories[v.key])
	end
	open(savedir * v.allele_freq_file, "w") do io
		writedlm(io, allele_frequencies[v.key])
	end
	open(savedir * v.switch_info_file, "w") do io
		CSV.write(io, switch_info[v.key])
	end
end

# ╔═╡ 344605b7-7735-49e2-b965-37dac463ab85
md"# Tests"

# ╔═╡ ab7f8ac2-720b-4c76-8db2-1118ab753ebf


# ╔═╡ 848b216e-3759-4e79-a6f3-63163c58a3fe


# ╔═╡ ec2f318a-082c-4213-8597-a85e491071e3


# ╔═╡ 2fda3321-ff5a-4134-8646-7da8b2666c6e


# ╔═╡ Cell order:
# ╠═10b7789c-e56f-4b71-8e38-312af11bb751
# ╠═19759f0e-ced3-11ed-2b40-8b01992454e2
# ╠═4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
# ╠═d8666482-d8a4-43f4-ad21-81f6e972ab15
# ╟─6cf21eab-ddb6-4c91-9e6a-d925d61702c2
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╟─c7909693-203b-49a8-8713-387dc287e4a2
# ╟─a37e9a1e-e06f-4402-abf3-f82e22b63fe3
# ╠═afe7bf04-d998-46fc-afed-8e2cd21361d2
# ╠═1110d7fa-9c96-4b40-ae60-2396a08d4db0
# ╠═f77e8aa9-7843-4917-a745-0fa3e6611939
# ╠═b7d9a8d0-476d-4d87-acd8-a06c7766198b
# ╠═420fceea-4103-4262-95f8-0fb4450da83f
# ╠═880c9804-d0e3-4f71-a34c-fcef988e3953
# ╟─51b03b84-3946-4079-a43a-7ae0050f2b93
# ╟─cb1a8bf0-776d-41d1-8a50-90c81292c5be
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═e8a8acd2-0aa0-4c36-b7be-8ef356fc6692
# ╠═50075364-ebd7-43b8-b1b0-3aa0375b01b4
# ╠═8172b776-f207-41fe-bd85-8fd2fada723f
# ╠═faf656e1-a817-4a30-a9e8-f86cb37ef35d
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═ab7f8ac2-720b-4c76-8db2-1118ab753ebf
# ╠═848b216e-3759-4e79-a6f3-63163c58a3fe
# ╠═ec2f318a-082c-4213-8597-a85e491071e3
# ╠═2fda3321-ff5a-4134-8646-7da8b2666c6e
