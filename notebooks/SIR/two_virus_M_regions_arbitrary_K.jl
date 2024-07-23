### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 24a8ea84-723e-11ee-2c9b-f39f058981c6
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using Distributions
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR
	using StatsBase
end

# ╔═╡ 7b9cdd9a-0d63-4d9a-b260-5c0d30561b2b
md"# Setup"

# ╔═╡ f6e50cb8-97c5-41a3-8c3c-2f39eecbe4e8
md"""
Two sliders: 
- `M` controls the number of regions;
- `c` controls the connectivity between regions, with $0 \leq c \leq 1/M$

The `C` matrix can be accesssed in `params.C`
"""

# ╔═╡ d349ee12-427b-4a50-9a9e-4db775a3b339
function log_range(low, high, length)
	@assert low > 0 && high > low
	exp.(range(log(low), log(high); length))
end

# ╔═╡ 88c1400c-e24e-4734-98e0-9d3d7f9a29a4
begin
	local _Ms = @bind M Slider([1,2,3,4,5,10,20], show_value=true, default=2)
	local c_values = vcat([0], log_range(1/M*1e-3, 1/M, 10))
	local _Cs = @bind c Slider(c_values, show_value=true, default = 0.)
	Ms = md"M = $(_Ms)"
	Cs = md"c = $(_Cs)"
end;

# ╔═╡ 404644b7-5a70-4fc4-a9ad-c0ad0f020ca3
Ms

# ╔═╡ b1b34106-ba8b-4e45-8a18-243275755edf
Cs

# ╔═╡ c6ab1a8a-fbf8-40f3-99d0-124190a5823e
params = let
	N = 2
	α = 3
	γ = .005
	PSS.Parameters(; N, M, α, γ, c)
end;

# ╔═╡ 058a199d-8e04-4321-87db-42129c20657c
@unpack α, γ, δ = params;

# ╔═╡ ddaf4d65-e698-4b80-83a2-885d3cc6eade
T = 15/params.γ # simulation time

# ╔═╡ 7383cc25-5b69-4c5b-9339-c833944bade9
md"Cross immunity matrix for the special region $i=1$."

# ╔═╡ 2d7ddbdc-7377-4700-8840-6e820b5486ed
ϕ_wt = 1

# ╔═╡ 6df8c075-e43a-4a8a-9ab3-cd1d1a2551d4
ϕs = let
	_slider = @bind ϕ_m Slider(.9:0.01:1.1, show_value=true, default=1)
	ϕs = md"""
	 $\phi_m$  = $(_slider) 
	"""
end

# ╔═╡ e113a279-e603-4be4-90bf-f105f53bca90
ξs = md" $\xi$ = $((ϕ_m - δ/α)/(ϕ_wt - δ/α))"

# ╔═╡ 11516882-5f70-4d7b-8838-e90df96a9c6b
ϕ_m_slider = (ϕs, ξs)

# ╔═╡ a1006bdd-8f85-4e43-9b1b-5317df16fe05
md"Cross immunity matrix for other regions $i>1$."

# ╔═╡ 7c390914-ac62-4108-9a7d-1c8bebb21efb
md"""
Setting up regions in a vector `R`: 
- `R[1]` has a special cross-immunity matrix
- other regions have `K0`
"""

# ╔═╡ 12c52ebc-c7f9-4cd6-b198-60a2e1af637e
md"# Simulation"

# ╔═╡ 5ff6d9d0-482b-4db0-b928-2e9f6f81211b
md"Simulating for time `T` to get the initial state, and then introducing the mutant in all regions."

# ╔═╡ fcce6223-9c61-4e0c-aaf8-0fccd33c26c5
md"## Quick note on indexing"

# ╔═╡ ddc1d57c-5092-417f-91d9-65bc30ce5139
md"""
The type `PartialSweepSIR.SIRState`, like the variables `state_init` and `state_final`, can be indexed in the following way: 
- `state[i, a, X]` where `X` is in `[:S, :I, :C, :R]` returns the size of compartment `X` in region `i` for virus `a`. `i` and `a` can also be ranges, *e.g.* `state[1, 1:2, :I]`.
- `state[:, a, X]` returns a vector of size `M` with the size of compartment `X` in all regions $i\in[1,M]$
- `state[i, :, X]`: same, but for all viruses
"""

# ╔═╡ 1dedef59-8645-494a-a92e-9b8b31a1b237
md"""
The type `PartialSweepSIR.Solution`, like the variable `sol`, can be used like this: 
- `sol(t)` returns the state at time `t` in the form of an `SIRState` object. The paragraph above shows how to index it.
- `sol[tvals, i, a, X]` returns the vector of the size of compartment `X` in region `i` for virus `a` and for all values in `tvals`.  
"""

# ╔═╡ c573178d-dc84-4c00-9438-86a25ca6d9ad
md"# Figures"

# ╔═╡ d07e4860-f7b4-427b-b7da-5ecc36731207
Ms

# ╔═╡ e5fe27b4-022e-4f00-a8bf-332f10679934
ϕ_m_slider

# ╔═╡ 4de36083-3848-4f89-b45d-47da15a30a10
Cs

# ╔═╡ eb61a918-14c2-4070-9544-0135adbe01f4
md"# Tests"

# ╔═╡ db0d5371-ab98-4181-9524-f564883d2f67
ϕs

# ╔═╡ f6855819-1064-4e57-9a8c-b17eedbe343d
Ms

# ╔═╡ 30723816-9fa4-4372-bfab-29ecfdf5d82d
Cs

# ╔═╡ 771b710d-ab4c-408f-92d5-165d6dcf3f57
Ks = begin
	ε = .1
	bs = 1 .- rand(Exponential(ε), M)
	fs = 1 .- rand(Exponential(ε), M)
	map(zip(bs, fs)) do (b, f)
		[1 b; f 1]
	end
end

# ╔═╡ 305e8857-8be5-4626-81f5-450d679249be
regions = let
	S0 = .4
	I0 = 1e-6
	C0 = 0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0, ϕ=ϕ_wt) # wt virus
	v2 = PSS.Virus(; S=S0, I=0, C=0., ϕ=ϕ_m) # mutant: initially not here
	
	R = [PSS.Region(; viruses=[v1,v2], K=Ks[m]) for m in 1:params.M]
	# r1 = PSS.Region(; viruses=[v1,v2], K=K_special)
	# R[1] = r1 # first region has the special cross-immunity
	# R
end

# ╔═╡ 9a140125-d656-4e2a-a856-b009dc8a09ce
state_init ,sol_init = let
	state = PSS.SIRState(; regions, parameters=params)
	sol = PSS.simulate(state, (0, T))
	PSS.set_infected(sol(T), 2, 1e-6), sol
end;

# ╔═╡ 97b68f57-e8c0-45d0-9915-1df738ea0153
sol = PSS.simulate(state_init, (0, T));

# ╔═╡ 2a64ac84-bc33-4ded-8a85-c8504b96e99a
let
	tvals = range(sol.tspan..., length=100)
	plot(PSS.frequency(sol, tvals, 2, 2) - PSS.frequency(sol, tvals, 1, 2))
end

# ╔═╡ 39e7c4a3-9cc2-4a6d-b726-dceba43007b0
f_av = let
	f_av = mean(fs)
	b_av = mean(bs)
	(1 - f_av)/(2 - f_av - b_av)
end

# ╔═╡ 603ca451-8d5c-4f60-91ff-8f1b198b95b0
p_freq = let
	# freq. plot
	tvals = range(sol.tspan..., length=100)
	f_R1 = PSS.frequency(sol, tvals, 1, 2)
	f_R2 = PSS.frequency(sol, tvals, 2, 2)
	f_R = PSS.frequency(sol, tvals, 2)

	lw = 2

	p = plot(
		legend = :bottomright,
		xlabel = "Time",
		ylabel = "",
		title = "Mutant frequency",
		ylim = (-0.03, 1.03),
	)
	
	# plot!(tvals, f_R1, label="Region 1", linewidth=lw)
	# plot!(tvals, f_R2, label="Other regions", linewidth=lw)
	for m in 1:M
		plot!(tvals, PSS.frequency(sol, tvals, m, 2), line=(:blue, lw/2, .1), label="")
	end
	plot!(tvals, f_R, label="Overall", linewidth=2*lw, color=1)

	hline!([f_av], label="", line=(:black))
	p
end

# ╔═╡ fb29611a-85aa-4c32-9771-f95dfc3f5f55
x_m_c0 = begin
	num = sum(zip(fs, bs)) do (f, b)
		(1-f)/(1-b*f)
	end
	den = sum(zip(fs, bs)) do (f, b)
		(2-f-b)/(1-b*f)
	end
	num/den
end

# ╔═╡ b5378128-efd4-409b-84e2-62f0b6689a8c
x_m_clarge = let
	Kav = mean(Ks)
	@info (Kav=Kav,)
	(1 - Kav[2,1]) / (2 - Kav[1,2] - Kav[2,1])
end

# ╔═╡ 56ead77b-db84-407e-b8eb-64d1639f0477
let
	# freq. plot
	tvals = range(sol.tspan..., length=100)
	f_R = PSS.frequency(sol, tvals, 2)

	lw = 2

	p = plot(
		legend = :bottomright,
		xlabel = "Time",
		ylabel = "",
		title = "Mutant frequency",
		ylim = (-0.03, 1.03),
	)
	for m in 1:M
		plot!(tvals, PSS.frequency(sol, tvals, m, 2), line=(:blue, 1, .3), label="")
	end
	# plot!(tvals, f_R, label="Overall", linewidth=1, color=1)

	# hline!([x_m_c0], label="eq. c=0", line=(:black))
	# hline!([x_R1], label="eq. R1", line=(:black, :dashdot))
	hline!([x_m_clarge], label="eq. c=1/M", line=(:black, :dash))
	p
end

# ╔═╡ 10f605f9-9a00-49eb-a673-7f485b8d615b
bs

# ╔═╡ 20dcc678-84f6-414c-a34f-36fc5a0b2f95
fs

# ╔═╡ dd8b3282-472f-4e55-a263-a6712bfd65f3
x_R1 = (1 - fs[1]) / (2 - fs[1] - bs[1])

# ╔═╡ bde27015-8d5f-4c99-afb4-8159c2f1cd6b
begin
	f_eff = fs[1]/M + (1-ε)*(M-1)/M
	b_eff = bs[1]/M + (1-ε)*(M-1)/M
	(1-f_eff)/(2 - b_eff - f_eff)
end

# ╔═╡ 19b3cd13-9d8d-4631-909a-7589b246ba7d
let
	tvals = range(sol.tspan..., length=100)
	i=1
	p = plot(tvals, sol[tvals, 1, 2, :S])
	plot!(tvals, sol[tvals, 2, 2, :S])
	plot!(tvals, sol[tvals, 3, 2, :S])

	Im = mean(m -> sol[tvals[end], m, 2, :I], 1:M)
	Kav = mean(Ks)
	I = γ/δ / (1-δ/α) * (Kav \ [1., 1.])
	# hline!([δ/α])
	# hline!([Sm])
	# hline!([mean(I)])
	p
end

# ╔═╡ c1e90581-7d9a-41de-9b1d-4c6129f56e16
let
	tvals = range(sol.tspan..., length=100)
	
	i = 2
	I = [mean(m -> sol[tvals[end], m, a, :I], 1:M) for a in 1:2]
	for a in 1:2
		@info "S_$(i)^$(a) = " sol[tvals[end], i, a, :S]
		@info "th" δ/α * sol[tvals[end], i, a, :I] / I[a]
	end
end

# ╔═╡ f1e17091-d5b6-4085-859e-a02928dbc6cb
let
	tvals = range(sol.tspan..., length=100)
	i=3
	Ii = [sol[tvals[end], i, a, :I] for a in 1:2]
	Ii_s = sum(Ii)
	Imat = [sol[tvals[end], i, a, :I] for i in 1:M, a in 1:2]
	Smat = [sol[tvals[end], i, a, :S] for i in 1:M, a in 1:2]
	Itot = sum(m -> sum(a -> sol[tvals[end], i, a, :I], 1:2), 1:M)
	
	D = map(1:2) do _
		γ/δ*(1 - M*δ/α * Ii_s/Itot)
	end

	@info Ii
	@info Ks[i] \ D
	@info Ii_s / Itot
	@info Ii

	Imat
	ν = Imat[:,2] ./ sum(Imat, dims=2) |> mean
	Imat[:,1]  - sum(Imat, dims=2)* (1-ν)

	Smat
	Smat .- M*δ/α * sum(Imat, dims=2) ./ Itot
end

# ╔═╡ f8f8a5e4-ace8-4dee-a10a-a72e8e35dc06
Ks

# ╔═╡ 0b663a1e-0c3c-4860-8de6-73f1aafcf0dc
md"# END TESTS"

# ╔═╡ 9d9c92da-1fb2-40cc-a7f7-b3025edefed7
pS = let
	g = :S
	
	tvals = collect(range(sol.tspan..., length=400))
	lw = .5
	
	X_wt = map(t -> mean(sol[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol[t, 1:M, 2, g]), tvals)
	
	p = plot(
		title="Susceptibles (all regions)", 
		xlabel="Time",
		legend=:topright,
	)
	plot!(tvals, X_wt, label="w.t.", linewidth = lw)
	plot!(tvals, X_m, label="Mutant", linewidth = lw)
	hline!([δ/α], line=(lw+2, :black, 0.3), label="")

	p
end

# ╔═╡ dd5255ef-3045-48f8-8b51-f97b0c006c75
pI = let
	g = :I
	tvals = range(sol.tspan..., length=400)
	lw = 2
	
	X_wt = map(t -> mean(sol[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol[t, 1:M, 2, g]), tvals)
	
	p = plot(
		title="Infectious (all regions)", 
		xlabel="Time",
		legend=:topright,
	)
	plot!(tvals, X_wt, label="w.t.", linewidth = lw)
	plot!(tvals, X_m, label="Mutant", linewidth = lw)

	p
end

# ╔═╡ Cell order:
# ╠═24a8ea84-723e-11ee-2c9b-f39f058981c6
# ╟─7b9cdd9a-0d63-4d9a-b260-5c0d30561b2b
# ╟─f6e50cb8-97c5-41a3-8c3c-2f39eecbe4e8
# ╠═d349ee12-427b-4a50-9a9e-4db775a3b339
# ╠═88c1400c-e24e-4734-98e0-9d3d7f9a29a4
# ╠═404644b7-5a70-4fc4-a9ad-c0ad0f020ca3
# ╠═b1b34106-ba8b-4e45-8a18-243275755edf
# ╠═c6ab1a8a-fbf8-40f3-99d0-124190a5823e
# ╠═058a199d-8e04-4321-87db-42129c20657c
# ╠═ddaf4d65-e698-4b80-83a2-885d3cc6eade
# ╟─7383cc25-5b69-4c5b-9339-c833944bade9
# ╠═2d7ddbdc-7377-4700-8840-6e820b5486ed
# ╟─e113a279-e603-4be4-90bf-f105f53bca90
# ╟─6df8c075-e43a-4a8a-9ab3-cd1d1a2551d4
# ╠═11516882-5f70-4d7b-8838-e90df96a9c6b
# ╟─a1006bdd-8f85-4e43-9b1b-5317df16fe05
# ╟─7c390914-ac62-4108-9a7d-1c8bebb21efb
# ╠═305e8857-8be5-4626-81f5-450d679249be
# ╟─12c52ebc-c7f9-4cd6-b198-60a2e1af637e
# ╟─5ff6d9d0-482b-4db0-b928-2e9f6f81211b
# ╠═9a140125-d656-4e2a-a856-b009dc8a09ce
# ╠═97b68f57-e8c0-45d0-9915-1df738ea0153
# ╟─fcce6223-9c61-4e0c-aaf8-0fccd33c26c5
# ╟─ddc1d57c-5092-417f-91d9-65bc30ce5139
# ╟─1dedef59-8645-494a-a92e-9b8b31a1b237
# ╟─c573178d-dc84-4c00-9438-86a25ca6d9ad
# ╠═d07e4860-f7b4-427b-b7da-5ecc36731207
# ╠═39e7c4a3-9cc2-4a6d-b726-dceba43007b0
# ╠═e5fe27b4-022e-4f00-a8bf-332f10679934
# ╠═4de36083-3848-4f89-b45d-47da15a30a10
# ╠═603ca451-8d5c-4f60-91ff-8f1b198b95b0
# ╟─eb61a918-14c2-4070-9544-0135adbe01f4
# ╠═fb29611a-85aa-4c32-9771-f95dfc3f5f55
# ╠═b5378128-efd4-409b-84e2-62f0b6689a8c
# ╠═10f605f9-9a00-49eb-a673-7f485b8d615b
# ╠═20dcc678-84f6-414c-a34f-36fc5a0b2f95
# ╠═dd8b3282-472f-4e55-a263-a6712bfd65f3
# ╠═bde27015-8d5f-4c99-afb4-8159c2f1cd6b
# ╠═db0d5371-ab98-4181-9524-f564883d2f67
# ╠═f6855819-1064-4e57-9a8c-b17eedbe343d
# ╠═30723816-9fa4-4372-bfab-29ecfdf5d82d
# ╠═56ead77b-db84-407e-b8eb-64d1639f0477
# ╠═2a64ac84-bc33-4ded-8a85-c8504b96e99a
# ╠═771b710d-ab4c-408f-92d5-165d6dcf3f57
# ╠═19b3cd13-9d8d-4631-909a-7589b246ba7d
# ╠═c1e90581-7d9a-41de-9b1d-4c6129f56e16
# ╠═f1e17091-d5b6-4085-859e-a02928dbc6cb
# ╠═f8f8a5e4-ace8-4dee-a10a-a72e8e35dc06
# ╟─0b663a1e-0c3c-4860-8de6-73f1aafcf0dc
# ╠═9d9c92da-1fb2-40cc-a7f7-b3025edefed7
# ╠═dd5255ef-3045-48f8-8b51-f97b0c006c75
