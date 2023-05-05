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
	using JSON3
	using Plots
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
βvals = [.1, .2 , .4, .6]

# ╔═╡ af5908c4-5183-4089-8647-069740dc13e9
[L*β^2 for β in βvals]

# ╔═╡ ef2209f0-6a67-4e28-a294-36a33c3faf6b
ρvals = let
	v = [20, 50, 100, 250, 1000]
	sort(1 ./ v)
end

# ╔═╡ 7ff23748-6c11-4fb9-ab62-69090ce37f59
s(β, α) = -α*log(1-β)/2

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
	[(ρ=ρ, β=β) => s(β, α)/2/ρ for β in βvals for ρ in ρvals], 
	(length(ρvals), length(βvals))
)

# ╔═╡ 0f0e6836-802a-4c9c-b274-54982a948030
evtime(ρ, β) = 75 / ρ / β^2

# ╔═╡ 09445e51-d5ee-4032-97e9-177f585244de
reshape(
	[evtime(ρ, β) for β in βvals for ρ in ρvals], 
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

# ╔═╡ 7e2650ad-406f-4553-b6a4-fea7e766f146
1/sqrt(L)

# ╔═╡ b33f4b77-ef9f-4c2c-be7d-0aade6087455
1 .- exp.(-ρvals/α)

# ╔═╡ 1927c41c-941e-4c2e-adec-23e2abe06777
ρvals[end]

# ╔═╡ bbb24ac7-d514-4343-8c07-5f017ce45e0b
function simulate(ρ, β)
	# setting parameters
	fitness_distribution = Dirac(s(β, α))
	T = evtime(ρ, β)
	Δt = round(Int, T/2_500) + 1 # 10_000 pts per simulation
	switchgen = round(Int, 1/ρ)

	# initial population
	H = zeros(Float64, L)
	ϕ = ExpiringFitness(; L, H, α)
	pop = Pop(ϕ; N, L, μ)
	cb = (
		# f1 = WF.frequencies,
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

	return cb_vals
end

# ╔═╡ 8531b7a7-83c9-44dd-b7c0-f69757d639ab


# ╔═╡ 2660ca6f-7554-4cf7-a322-a533418bd7b7
data = let
	dat = Dict()
	for ρ in ρvals, β in βvals
		@info ρ, β
		cb = simulate(ρ, β)
		df = DataFrame()
		for x in cb
			push!(df, x)
		end
		dat[ρ, β] = df
	end
	dat
end

# ╔═╡ 3164f91e-c112-4947-8670-d0b65d7a9667
filenames = Dict(
	"rho$(round(ρ; sigdigits=2))_beta$(round(β; sigdigits=2))_dat.csv" => (ρ=ρ, β=β, L=L, α=α)
	for (ρ, β) in keys(data)
)

# ╔═╡ 8172b776-f207-41fe-bd85-8fd2fada723f
open(savedir * "/files.json", "w") do io
	JSON3.pretty(io, filenames)
end

# ╔═╡ faf656e1-a817-4a30-a9e8-f86cb37ef35d
for (file, (ρ, β)) in filenames
	open(savedir * file, "w") do io
		CSV.write(io, data[ρ, β])
	end
end

# ╔═╡ 344605b7-7735-49e2-b965-37dac463ab85
md"# Tests"

# ╔═╡ 7f55239b-c970-4649-b920-1e45672f1913
cb = simulate(ρvals[1], βvals[end])

# ╔═╡ 65fb5f1d-8816-439d-9994-8b58c6f4f255
f1 = vcat(map(x -> x.f1', cb)...);

# ╔═╡ 85acf922-7f16-4115-a10b-9eeab7d76326
# ╠═╡ disabled = true
#=╠═╡
βvals = [0.8]
  ╠═╡ =#

# ╔═╡ d2f2fa3a-7fab-43fa-9f7a-56fa1f908d5b
1/ρvals[1]

# ╔═╡ 2907f8ae-ad60-45f2-a58e-453c7e87af8e
let
	plot(f1[:, 2:2:end], label="")
	plot!(xlim = 000 .+ (0,1000))
	# hline!([β1], color=:black)
end

# ╔═╡ 942e5565-0f49-4ffc-a70c-c94492e7f15a
log(30)

# ╔═╡ 6d51bff5-301c-42c9-9e80-c7adc8efd9d8
let
	x = map(v -> v.varpos_5, cb)
	plot(x*L)
	hline!([βvals[end]^-2])
end

# ╔═╡ 7ce0beaa-9227-4d99-af57-12c719cf0be5
let
	x = map(v -> v.varpos_strict, cb)
	plot(x*L)
end

# ╔═╡ 29b634a9-3a2e-4e60-adda-bfdfcd6ac24c
let
	x = map(v -> v.n_genotypes, cb)
	plot(x*L)
end

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
# ╠═7ff23748-6c11-4fb9-ab62-69090ce37f59
# ╟─51b03b84-3946-4079-a43a-7ae0050f2b93
# ╠═36a7ed0a-7a51-42c6-9f98-050f56be59ab
# ╠═0f0e6836-802a-4c9c-b274-54982a948030
# ╠═09445e51-d5ee-4032-97e9-177f585244de
# ╟─cb1a8bf0-776d-41d1-8a50-90c81292c5be
# ╠═7e2650ad-406f-4553-b6a4-fea7e766f146
# ╠═b33f4b77-ef9f-4c2c-be7d-0aade6087455
# ╠═1927c41c-941e-4c2e-adec-23e2abe06777
# ╠═bbb24ac7-d514-4343-8c07-5f017ce45e0b
# ╠═8531b7a7-83c9-44dd-b7c0-f69757d639ab
# ╠═2660ca6f-7554-4cf7-a322-a533418bd7b7
# ╠═3164f91e-c112-4947-8670-d0b65d7a9667
# ╠═8172b776-f207-41fe-bd85-8fd2fada723f
# ╠═faf656e1-a817-4a30-a9e8-f86cb37ef35d
# ╠═344605b7-7735-49e2-b965-37dac463ab85
# ╠═7f55239b-c970-4649-b920-1e45672f1913
# ╠═65fb5f1d-8816-439d-9994-8b58c6f4f255
# ╠═85acf922-7f16-4115-a10b-9eeab7d76326
# ╠═d2f2fa3a-7fab-43fa-9f7a-56fa1f908d5b
# ╠═2907f8ae-ad60-45f2-a58e-453c7e87af8e
# ╠═942e5565-0f49-4ffc-a70c-c94492e7f15a
# ╠═6d51bff5-301c-42c9-9e80-c7adc8efd9d8
# ╠═7ce0beaa-9227-4d99-af57-12c719cf0be5
# ╠═29b634a9-3a2e-4e60-adda-bfdfcd6ac24c
# ╠═ab7f8ac2-720b-4c76-8db2-1118ab753ebf
# ╠═848b216e-3759-4e79-a6f3-63163c58a3fe
# ╠═ec2f318a-082c-4213-8597-a85e491071e3
# ╠═2fda3321-ff5a-4134-8646-7da8b2666c6e
