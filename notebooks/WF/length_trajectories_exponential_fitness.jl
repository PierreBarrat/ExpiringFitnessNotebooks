### A Pluto.jl notebook ###
# v0.19.22

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
	L = 100
	μ = 0
	α = 0.1
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

# ╔═╡ c5d0b72c-a90a-4e05-8416-a333b85c23c8
βvals = [.15, .2 , .4, .6]

# ╔═╡ af5908c4-5183-4089-8647-069740dc13e9
[L*β^2 for β in βvals]

# ╔═╡ ef2209f0-6a67-4e28-a294-36a33c3faf6b
ρvals = let
	v = [10, 20, 50, 100, 250]
	sort(1 ./ v)
end

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

# ╔═╡ 43062de6-5809-4cf8-8843-c624860abdec
function get_Ne(ρ, β, α)
	s = get_s(β, α)
	β2 = 2*s^2 / (s+α) / (2*s+α)
	return 1/ρ/β2
end

# ╔═╡ 4a491e08-4ac9-462d-9677-78079c875aea
get_Δt(ρ, β, α) = round(Int, get_Ne(ρ, β, α) / 500) + 1

# ╔═╡ 51b03b84-3946-4079-a43a-7ae0050f2b93
md"""
An important ratio is $s/\rho$: 
- if $s \gg \rho$, the sweeps happen one after the other;
- if $s/\rho \ll 1$, they overlap a lot. 

As a function of $\beta$, we have

$$s/\rho = -α\rho^{-1}\cdot\log(1-\beta).$$
"""

# ╔═╡ 36a7ed0a-7a51-42c6-9f98-050f56be59ab
reshape(
	[(ρ=ρ, β=β) => get_s(β, α)/2/ρ for β in βvals for ρ in ρvals], 
	(length(ρvals), length(βvals))
)

# ╔═╡ 09445e51-d5ee-4032-97e9-177f585244de
reshape(
	[get_Ne(ρ, β, α) for β in βvals for ρ in ρvals],
	(length(ρvals), length(βvals))
)

# ╔═╡ cb1a8bf0-776d-41d1-8a50-90c81292c5be
md"""

It could be interesting to see what happens in the regime $s/ρ < 1$ and $Lβ^2 \gg 1$, meaning that we have a large supply of mutations according but the sweeps overlap. 

In terms of $β$ and $ρ$ this regime is defined by 

$$\beta < 1 - e^{-\rho/\alpha} \;\;\text{and}\;\; \beta \gg L^{-1/2}$$

This is of course only possible if 

$$1 - e^{-\rho/\alpha} \gg L^{-1/2}$$
"""

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(ρ, β; α=0., L=0, N=0, μ=0.)
	# setting parameters
	s = get_s(β, α)
	fitness_distribution = Exponential(s/2)

	T = max(500/ρ, 20 * get_Ne(ρ, β, α))
	Δt = get_Δt(ρ, β, α)
	
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
		change_field_time = :periodic,
		max_freq=1/N, 
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

# ╔═╡ 2660ca6f-7554-4cf7-a322-a533418bd7b7
trajectories, diversity = let
	tjs = Dict() 
	div = Dict()
	for ρ in ρvals, β in βvals
		@info ρ, β
		@time T, d = simulate(ρ, β; L, N, α, μ)
		tjs[ρ, β] = T
		
		df = DataFrame()
		for x in d
			push!(df, x)
		end
		div[ρ, β] = df
	end
	tjs, div
end

# ╔═╡ 3164f91e-c112-4947-8670-d0b65d7a9667
filenames = map(collect(keys(trajectories))) do (ρ, β)
	rβ = round(β, sigdigits=2)
	rρ = round(ρ, sigdigits=2)
	"rho$(rρ)_beta$(rβ)" => (
		ρ=ρ, β=β, L=L, α=α, N=N, 
		Ne=get_Ne(ρ, β, α), Δt=get_Δt(ρ, β, α), s = get_s(β, α),
		trajectory_file = "trajectory_rho$(rρ)_beta$(rβ).csv",
		diversity_file = "diversity_rho$(rρ)_beta$(rβ).csv",
	)
end |> Dict

# ╔═╡ 8172b776-f207-41fe-bd85-8fd2fada723f
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, filenames)
end

# ╔═╡ faf656e1-a817-4a30-a9e8-f86cb37ef35d
for v in values(filenames)
	β = v.β
	ρ = v.ρ
	open(savedir * v.diversity_file, "w") do io
		CSV.write(io, diversity[ρ, β])
	end
	open(savedir * v.trajectory_file, "w") do io
		CSV.write(io, trajectories[ρ, β])
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
# ╠═c5d0b72c-a90a-4e05-8416-a333b85c23c8
# ╠═af5908c4-5183-4089-8647-069740dc13e9
# ╠═ef2209f0-6a67-4e28-a294-36a33c3faf6b
# ╟─a75a63dd-8654-4eff-a95a-15bac2296c47
# ╠═7ff23748-6c11-4fb9-ab62-69090ce37f59
# ╠═43062de6-5809-4cf8-8843-c624860abdec
# ╠═4a491e08-4ac9-462d-9677-78079c875aea
# ╟─51b03b84-3946-4079-a43a-7ae0050f2b93
# ╠═36a7ed0a-7a51-42c6-9f98-050f56be59ab
# ╠═09445e51-d5ee-4032-97e9-177f585244de
# ╟─cb1a8bf0-776d-41d1-8a50-90c81292c5be
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═2660ca6f-7554-4cf7-a322-a533418bd7b7
# ╠═3164f91e-c112-4947-8670-d0b65d7a9667
# ╠═8172b776-f207-41fe-bd85-8fd2fada723f
# ╠═faf656e1-a817-4a30-a9e8-f86cb37ef35d
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═ab7f8ac2-720b-4c76-8db2-1118ab753ebf
# ╠═848b216e-3759-4e79-a6f3-63163c58a3fe
# ╠═ec2f318a-082c-4213-8597-a85e491071e3
# ╠═2fda3321-ff5a-4134-8646-7da8b2666c6e
