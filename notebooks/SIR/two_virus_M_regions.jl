### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ 4982e04e-6501-11ed-132e-ad3e395233da
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR
	using StatsBase
end

# ╔═╡ db663438-8523-47ac-a783-88ff59a7ab1d
md"# Setup"

# ╔═╡ c6c407d9-fea7-4e79-97ff-6acb88f4356d
md"""
Two sliders: 
- `M` controls the number of regions;
- `c` controls the connectivity between regions, with $0 \leq c \leq 1/M$

The `C` matrix can be accesssed in `params.C`
"""

# ╔═╡ c79fe879-a32a-4391-8df1-5b9d9eeb75f5
function log_range(low, high, length)
	@assert low > 0 && high > low
	exp.(range(log(low), log(high); length))
end

# ╔═╡ d42e08dd-198d-4b85-93fd-d9e522279b56
begin
	local _Ms = @bind M Slider([1,2,3,4,5,10,20], show_value=true, default=2)
	Ms = md"M = $(_Ms)"
end

# ╔═╡ 6459a142-ece5-4ad4-946e-708b684993c9
begin
	local c_values = vcat([0], log_range(1/M*1e-3, 1/M, 10))
	local _Cs = @bind c Slider(c_values, show_value=true, default = c_values[end])
	Cs = md"c = $(_Cs)"
end

# ╔═╡ f431ba40-9742-4898-92f0-33c3b8e12a1c
params = let
	N = 2
	α = 3
	γ = .005
	PSS.Parameters(; N, M, α, γ, c)
end;

# ╔═╡ e38f872f-f9ec-45cf-953c-42a18aa0bab3
@unpack α, γ, δ = params;

# ╔═╡ c2a159f4-cd25-4525-a270-2b3065a0dc0e
T = 5/params.γ # simulation time

# ╔═╡ a7fc518a-f6d3-4e93-bac5-02d946c768f8
md"Cross immunity matrix for the special region $i=1$."

# ╔═╡ c4da46f5-8178-4ae5-b865-db025a365da0
K_special = begin
	b_special = .7
	f_special = .8
	[1 b_special; f_special 1]
end

# ╔═╡ 81cbea6e-9a97-4e7e-9528-f0830c7ae58f
md"Cross immunity matrix for other regions $i>1$."

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
K0 = begin
	x = .97
	b = x
	f = x
	[1 b; f 1]
end

# ╔═╡ b3f03e99-9550-4b01-bda1-b7f9f6617c7d
Ms

# ╔═╡ b14f025b-e95f-42fb-a321-4d914f6072ba
Cs

# ╔═╡ b7dd719a-118e-4809-84da-609735534835
md"""
Setting up regions in a vector `R`: 
- `R[1]` has a special cross-immunity matrix
- other regions have `K0`
"""

# ╔═╡ a62f2142-7d76-4a15-9566-54115aa98e08
regions = let
	S0 = .4
	I0 = 1e-6
	C0 = 0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0) # wt virus
	v2 = PSS.Virus(; S=S0, I=0, C=0.) # mutant: initially not here
	
	R = [PSS.Region(; viruses=[v1,v2], K=K0) for m in 1:params.M]
	r1 = PSS.Region(; viruses=[v1,v2], K=K_special)
	R[1] = r1 # first region has the special cross-immunity
	R
end

# ╔═╡ 51c77696-d4ae-466e-a19a-cf9278c7d3a5
md"# Simulation"

# ╔═╡ 559e0978-3bb1-407f-b27b-4e31dc9df04e
md"Simulating for time `T` to get the initial state, and then introducing the mutant in all regions."

# ╔═╡ fd0d5424-074d-4fee-a41c-a80e61591348
state_init ,sol_init = let
	state = PSS.SIRState(; regions, parameters=params)
	sol = PSS.simulate(state, (0, T))
	PSS.set_infected(sol(T), 2, 1e-6), sol
end;

# ╔═╡ 62c7a1db-42d0-463f-9f3d-954d36f1b1d2
sol = PSS.simulate(state_init, (0, T));

# ╔═╡ 64c6f22b-9717-46c4-93cf-e1eafd107931
md"## Quick note on indexing"

# ╔═╡ 1d3084e3-a18e-4850-b424-a7a32f4780b3
md"""
The type `PartialSweepSIR.SIRState`, like the variables `state_init` and `state_final`, can be indexed in the following way: 
- `state[i, a, X]` where `X` is in `[:S, :I, :C, :R]` returns the size of compartment `X` in region `i` for virus `a`. `i` and `a` can also be ranges, *e.g.* `state[1, 1:2, :I]`.
- `state[:, a, X]` returns a vector of size `M` with the size of compartment `X` in all regions $i\in[1,M]$
- `state[i, :, X]`: same, but for all viruses
"""

# ╔═╡ 7aa5f567-8af0-4b40-8669-f959da124826
md"""
The type `PartialSweepSIR.Solution`, like the variable `sol`, can be used like this: 
- `sol(t)` returns the state at time `t` in the form of an `SIRState` object. The paragraph above shows how to index it.
- `sol[tvals, i, a, X]` returns the vector of the size of compartment `X` in region `i` for virus `a` and for all values in `tvals`.  
"""

# ╔═╡ bb4a3499-958f-4b34-b69c-de2874ae8546
md"# Figures"

# ╔═╡ 88b4146d-d8bb-4793-bec2-76a4cbc4e320
Ms

# ╔═╡ d17bc99f-f5cb-49c5-bfe1-e9219b64d19d
Cs

# ╔═╡ 8e17f448-7e6c-493c-b099-1f639a56a575
f_av = let
	f_av = (M-1)/M * f + f_special/M
	b_av = (M-1)/M * b + b_special/M
	(1 - f_av)/(2 - f_av - b_av)
end

# ╔═╡ 95efa1bd-8fe2-41f1-ba13-839a670cd86e
let
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
	
	plot!(tvals, f_R1, label="Region 1", linewidth=lw)
	plot!(tvals, f_R2, label="Other regions", linewidth=lw)
	plot!(tvals, f_R, label="Overall", linewidth=lw)
	# hline!(
	# 	[(1-f_special)/(2-b_special-f_special)]; 
	# 	label="", line=(:black, lw+2, 0.3)
	# )
	hline!([f_av], label="", line=(:black))
	p
end

# ╔═╡ 87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
p_freq = let
	# freq. plot
	tvals = range(sol.tspan..., length=100)
	# f_R1 = PSS.frequency(sol, tvals, 1, 2)
	# f_R2 = PSS.frequency(sol, tvals, 2, 2)
	f_R = PSS.frequency(sol, tvals, 2)

	lw = 2

	p = plot(
		legend = :bottomright,
		xlabel = "Time",
		ylabel = "",
		title = "Mutant frequency",
		ylim = (-0.03, 1.03),
		frame=:box
	)
	
	# plot!(tvals, f_R1, label="Region 1", linewidth=lw)
	# plot!(tvals, f_R2, label="Other regions", linewidth=lw)
	plot!(tvals, f_R, label="", linewidth=lw)
	# hline!(
	# 	[(1-f_special)/(2-b_special-f_special)]; 
	# 	label="", line=(:black, lw+2, 0.3)
	# )
	hline!([f_av], label="", line=(:black))

	savefig("/home/pierrebc/Documents/BaleLabo/Slides/LJP_2023/mfreq_M10.png")
	
	p
end

# ╔═╡ 377398c1-fa28-476c-9328-5ea26fadbb7e
Ms

# ╔═╡ 399b3907-a990-409d-8147-de87f52b5c37
Cs

# ╔═╡ add14358-6d60-47a5-876f-8e0651781c57
pS = let
	g = :S
	
	tvals = collect(range(sol.tspan..., length=400))
	lw = 2
	
	X_wt = map(t -> mean(sol[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol[t, 1:M, 2, g]), tvals)
	
	p = plot(
		title="Susceptibles (all groups)", 
		xlabel="Time",
		legend=:topright,
		frame=:box,
	)
	plot!(tvals, X_wt, label="w.t.", linewidth = lw)
	plot!(tvals, X_m, label="Mutant", linewidth = lw)
	hline!([δ/α], line=(lw+2, :black, 0.3), label="")

	savefig("/home/pierrebc/Documents/BaleLabo/Slides/LJP_2023/S_M10.png")

	p
end

# ╔═╡ 2e5f1f19-d5f2-4e54-8337-49db178b13c9
pI = let
	g = :I
	tvals = range(sol.tspan..., length=400)
	lw = 2
	
	X_wt = map(t -> mean(sol[t, 1:M, 1, g]), tvals)
	X_m = map(t -> mean(sol[t, 1:M, 2, g]), tvals)
	
	p = plot(
		title="Infectious (all groups)", 
		xlabel="Time",
		legend=:topright,
		frame=:box,
	)
	plot!(tvals, X_wt, label="w.t.", linewidth = lw)
	plot!(tvals, X_m, label="Mutant", linewidth = lw)

	savefig("/home/pierrebc/Documents/BaleLabo/Slides/LJP_2023/I_M10.png")
	
	p
end

# ╔═╡ 2a1c18b3-587c-47d6-9def-1794fa2e9889
let
	plot(pS, pI, layout=(1,2))
	savefig("/home/pierrebc/Documents/BaleLabo/Slides/LJP_2023/SI_M10.png")
end

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╟─db663438-8523-47ac-a783-88ff59a7ab1d
# ╟─c6c407d9-fea7-4e79-97ff-6acb88f4356d
# ╠═c79fe879-a32a-4391-8df1-5b9d9eeb75f5
# ╠═d42e08dd-198d-4b85-93fd-d9e522279b56
# ╠═6459a142-ece5-4ad4-946e-708b684993c9
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═e38f872f-f9ec-45cf-953c-42a18aa0bab3
# ╠═c2a159f4-cd25-4525-a270-2b3065a0dc0e
# ╟─a7fc518a-f6d3-4e93-bac5-02d946c768f8
# ╠═c4da46f5-8178-4ae5-b865-db025a365da0
# ╟─81cbea6e-9a97-4e7e-9528-f0830c7ae58f
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╠═b3f03e99-9550-4b01-bda1-b7f9f6617c7d
# ╠═b14f025b-e95f-42fb-a321-4d914f6072ba
# ╠═95efa1bd-8fe2-41f1-ba13-839a670cd86e
# ╟─b7dd719a-118e-4809-84da-609735534835
# ╠═a62f2142-7d76-4a15-9566-54115aa98e08
# ╟─51c77696-d4ae-466e-a19a-cf9278c7d3a5
# ╟─559e0978-3bb1-407f-b27b-4e31dc9df04e
# ╠═fd0d5424-074d-4fee-a41c-a80e61591348
# ╠═62c7a1db-42d0-463f-9f3d-954d36f1b1d2
# ╟─64c6f22b-9717-46c4-93cf-e1eafd107931
# ╟─1d3084e3-a18e-4850-b424-a7a32f4780b3
# ╟─7aa5f567-8af0-4b40-8669-f959da124826
# ╟─bb4a3499-958f-4b34-b69c-de2874ae8546
# ╠═88b4146d-d8bb-4793-bec2-76a4cbc4e320
# ╠═d17bc99f-f5cb-49c5-bfe1-e9219b64d19d
# ╠═8e17f448-7e6c-493c-b099-1f639a56a575
# ╠═87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
# ╠═377398c1-fa28-476c-9328-5ea26fadbb7e
# ╠═399b3907-a990-409d-8147-de87f52b5c37
# ╠═add14358-6d60-47a5-876f-8e0651781c57
# ╠═2e5f1f19-d5f2-4e54-8337-49db178b13c9
# ╠═2a1c18b3-587c-47d6-9def-1794fa2e9889
