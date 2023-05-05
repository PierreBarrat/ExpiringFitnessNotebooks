### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 10b7789c-e56f-4b71-8e38-312af11bb751
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
end

# ╔═╡ 19759f0e-ced3-11ed-2b40-8b01992454e2
notebook_name = replace(@__FILE__, r"\.jl#==#.*" => "") |> basename

# ╔═╡ d8666482-d8a4-43f4-ad21-81f6e972ab15
savedir = "data_" * notebook_name * "/"

# ╔═╡ 4a8fa33e-f996-4cf9-bd95-05de2bbd11f5
mkpath(savedir)

# ╔═╡ 0ed067a4-3780-432b-bd54-9d52a22f2bbf
begin
	N = 1_000_000
	L = 150
	μ = 0
	# α = 0.1
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
	Ne_vals = [100, 500, 2000]
	βvals = [.2 , .4, .6]
	αvals = [.01, .025, .1]
	Δtvals = [1, 3, 10, 30]
	md"**Variable parameters**"
end

# ╔═╡ 1110d7fa-9c96-4b40-ae60-2396a08d4db0
parameters = map(Iterators.product(Δtvals, αvals, βvals, Ne_vals)) do (Δt, α, β, Ne)
	(Ne=Ne, β=β, α=α, Δt=Δt)
end |> (x -> vcat(x...));

# ╔═╡ a75a63dd-8654-4eff-a95a-15bac2296c47
md"""
Assume fitness effects are exponentially distributed: 

$$P(s) \propto e^{-s / s_0}.$$

One can show that $\beta$ then follows a Beta distribution $B(1, \alpha/s_0)$: 

$$P(\beta) \propto (1-\beta)^{\alpha/s_0 - 1},$$

with the following moments: 

$$\langle \beta \rangle = \frac{s_0}{\alpha + s_0} \;\;\text{and}\;\; \langle\beta^2\rangle = \frac{s_0}{s_0 + \alpha} \frac{2s_0}{2 s_0 + \alpha}.$$

If we want to fix $\beta^2$, we want

$$\begin{align}
2(1 - \langle \beta^2\rangle) s_0^2 -\alpha\langle\beta^2\rangle s_0  - \alpha^2\langle\beta^2\rangle = 0, \;\;\text{etc...}\\
\end{align}$$


But it is easier to fix $\langle \beta \rangle$, one should thus pick $s_0$ such that

$$s_0 = \alpha\frac{\langle\beta\rangle}{1 - \langle\beta\rangle}.$$

Since the Wright-Fisher code uses $\{-1,1\}$ states, we divide $s$ obtained by two before running the simulation.
"""

# ╔═╡ 7ff23748-6c11-4fb9-ab62-69090ce37f59
get_s(β, α)  = α * β / (1-β)

# ╔═╡ a387a506-d0e6-4c11-820d-a6745337705f
function get_ρ(β, α, Ne)
	s = get_s(β, α)
	β2 = 2*s^2 / (s+α) / (2*s+α)
	return 1/Ne/β2
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
function simulate(Ne, β, α; Δt = 1, L=0, N=0, μ=0.)
	# setting parameters
	s = get_s(β, α)
	fitness_distribution = Exponential(s/2)

	ρ = get_ρ(β, α, Ne)
	T = max(500/ρ, 50 * Ne)
	
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

	return trajectories, diversity
end

# ╔═╡ e8a8acd2-0aa0-4c36-b7be-8ef356fc6692
trajectories, diversity = let
	tjs = Dict() 
	div = Dict()
	for (i, p) in enumerate(parameters)
		@info p i/length(parameters)
		@time T, d = simulate(p.Ne, p.β, p.α; L, N, μ, Δt=p.Δt)
		tjs[p] = T
		
		df = DataFrame()
		for x in d
			push!(df, x)
		end
		div[p] = df
	end
	tjs, div
end;

# ╔═╡ 50075364-ebd7-43b8-b1b0-3aa0375b01b4
filenames = map(enumerate(collect(keys(trajectories)))) do (idx, p)
	idx => (
		Ne = p.Ne,
		β = p.β,
		α = p.α,
		Δt = p.Δt,
		s = get_s(p.β, p.α),
		ρ = get_ρ(p.β, p.α, p.Ne),
		L=L, 
		N=N,
		trajectory_file = "trajectory_$(idx).csv",
		diversity_file = "diversity_$(idx).csv",
		idx = idx,
		key = p,
	)
end

# ╔═╡ 8172b776-f207-41fe-bd85-8fd2fada723f
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, filenames)
end

# ╔═╡ faf656e1-a817-4a30-a9e8-f86cb37ef35d
for (i, v) in filenames
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[v.key])
	end
	open(savedir * v.trajectory_file, "w") do io
		CSV.write(io, trajectories[v.key])
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
# ╠═0ed067a4-3780-432b-bd54-9d52a22f2bbf
# ╟─c7909693-203b-49a8-8713-387dc287e4a2
# ╟─a37e9a1e-e06f-4402-abf3-f82e22b63fe3
# ╠═afe7bf04-d998-46fc-afed-8e2cd21361d2
# ╠═1110d7fa-9c96-4b40-ae60-2396a08d4db0
# ╟─a75a63dd-8654-4eff-a95a-15bac2296c47
# ╠═7ff23748-6c11-4fb9-ab62-69090ce37f59
# ╠═a387a506-d0e6-4c11-820d-a6745337705f
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
