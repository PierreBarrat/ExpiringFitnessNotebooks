### A Pluto.jl notebook ###
# v0.19.25

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

# ╔═╡ 287520e8-97a3-485b-9634-fd2db7a2f4e6
begin
	local _Ms = @bind M Slider([1,2,3,4,5,10,20], show_value=true, default=2)
	Ms = md"M = $(_Ms)"
end;

# ╔═╡ 3074735c-5073-4a7d-8cc9-d983cc0268f2
begin
		local c_values = vcat([0], log_range(1/M*1e-3, 1/M, 10))
		local _Cs = @bind c Slider(c_values, show_value=true, default = 1/M)
		Cs = md"c = $(_Cs)"
end

# ╔═╡ b3f03e99-9550-4b01-bda1-b7f9f6617c7d
Ms

# ╔═╡ b14f025b-e95f-42fb-a321-4d914f6072ba
Cs

# ╔═╡ f431ba40-9742-4898-92f0-33c3b8e12a1c
params = let
	N = 2
	α = 3
	γ = .005
	PSS.Parameters(; N, M, α, γ, c)
end;

# ╔═╡ e38f872f-f9ec-45cf-953c-42a18aa0bab3
@unpack α, γ, δ = params;

# ╔═╡ 644eed93-fdcc-479d-823d-076d40b517d5
params.C

# ╔═╡ c2a159f4-cd25-4525-a270-2b3065a0dc0e
T = 50/params.γ # simulation time

# ╔═╡ bd7c5323-1763-411a-9d94-ba25d77ea8d9
md"### Cross-immunity"

# ╔═╡ a7fc518a-f6d3-4e93-bac5-02d946c768f8
md"Cross immunity matrix for the special region $i=1$."

# ╔═╡ c4da46f5-8178-4ae5-b865-db025a365da0
K_special = begin
	b_special = .8
	f_special = .6
	[1 b_special; f_special 1]
end

# ╔═╡ 81cbea6e-9a97-4e7e-9528-f0830c7ae58f
md"Cross immunity matrix for other regions $i>1$."

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
K0 = begin
	b = .95
	f = .95
	[1 b; f 1]
end

# ╔═╡ 79bc278a-2e25-4ee4-b303-2f1f09cc3eac
md"### Fitness"

# ╔═╡ 8ff3efd3-0771-4e84-b700-4177f90f72cb
ϕ_wt = 1

# ╔═╡ c26c16b7-af52-46ec-9013-b7c353876262
ϕs = let
	_slider = @bind ϕ_m Slider(.9:0.01:1.1, show_value=true, default=1)
	ϕs = md"""
	 $\phi_m$  = $(_slider) 
	"""
end

# ╔═╡ f8025972-761d-4354-b688-35edaa639a80
ξs = md" $\xi$ = $((ϕ_m - δ/α)/(ϕ_wt - δ/α))"

# ╔═╡ 8e392fa9-4f6d-4b32-9c39-8e5c5000beb3
ϕ_m_slider = (ϕs, ξs)

# ╔═╡ abab9320-5e1f-475b-aa65-09c79f4cc284
ξ = (ϕ_m - δ/α)/(ϕ_wt - δ/α)

# ╔═╡ 2bf1ad78-0389-445b-8334-c6f9bea86dcf
md"""
The co-existence of the two strains depends on 
- the cross immunity $\epsilon_b$ and $\epsilon_f$;
- the ratio of growth rates $\xi=(\alpha\phi^m - \delta)/(\alpha\phi^{wt}-\delta)$. If $\xi<1$, the mutant strain has a fitness detriment.  
"""

# ╔═╡ 8b265670-6152-409a-809a-f68c45203196
md"""
### Co-existence condition: $c=0$

We just compute the average of the local equilibriums. Assuming $\xi < 1$, we have three cases: 

-  $\xi > \epsilon_f$

The mutant invades all regions, and we have the following equilibrium: 

$$\begin{align}
I^{wt} &= \frac{\gamma}{\delta}\frac{\alpha\phi^{wt} - \delta}{\alpha\phi^{wt}}\left(\frac{M-1}{M}\frac{1-\epsilon_b\xi}{1-\epsilon_b\epsilon_f} + \frac{1}{M}\frac{1-\epsilon_b^1\xi}{1-\epsilon_b^1\epsilon_f^1}\right)\\
I^{m} &= \frac{\gamma}{\delta}\frac{\alpha\phi^{m} - \delta}{\alpha\phi^{m}}\left(\frac{M-1}{M}\frac{1-\epsilon_f\xi^{-1}}{1-\epsilon_b\epsilon_f} + \frac{1}{M}\frac{1-\epsilon_f^1\xi^{-1}}{1-\epsilon_b^1\epsilon_f^1}\right)\\
\end{align}$$

-  $\epsilon_f^1 < \xi < \epsilon_f < 1$

In this case, the mutant can only invade the special region, and is absent from the others: 

$$\begin{align}
I^{wt} &= \frac{\gamma}{\delta}\frac{\alpha\phi^{wt} - \delta}{\alpha\phi^{wt}}\left(\frac{M-1}{M} + \frac{1}{M}\frac{1-\epsilon_b^1\xi}{1-\epsilon_b^1\epsilon_f^1}\right)\\
I^{m} &= \frac{1}{M}\frac{\gamma}{\delta}\frac{\alpha\phi^{m} - \delta}{\alpha\phi^{m}}\frac{1-\epsilon_f^1\xi^{-1}}{1-\epsilon_b^1\epsilon_f^1}\\
\end{align}$$



-  $\xi < \epsilon_f^1.$

In this case, the mutant cannot even invade the special region and $I_m = 0$.

"""

# ╔═╡ 9a14582a-1d87-49d6-8377-6df2a309da2d
begin
	I_eq_wt_c0 = γ/δ * (α*ϕ_wt - δ)/(α*ϕ_wt) * (
		(M-1)/M * min(1, max(0, (1-b*ξ)/(1-b*f))) +
		1/M * min(1, max(0, (1-b_special*ξ)/(1-b_special*f_special)))
	)
	I_eq_m_c0 = γ/δ * (α*ϕ_m - δ)/(α*ϕ_m) * (
		(M-1)/M * min(1, max(0, (1-f/ξ)/(1-b*f))) +
		1/M * min(1, max(0, (1-f_special/ξ)/(1-b_special*f_special)))
	)
	freq_eq_c0 = I_eq_m_c0 / (I_eq_m_c0 + I_eq_wt_c0)
end

# ╔═╡ 4643009b-aef8-469d-b25f-35df8b0593a9
md"""
### Co-existence condition: large $c$

We then have one effective region with 

$$\langle\epsilon_f\rangle = \frac{1}{M}\sum_i\epsilon_f^i$$ 

(and resp. for $\epsilon_b$). Assuming $\xi < 1$, the co-existence condition is 

$$\xi > \langle \epsilon_f \rangle.$$

We have two cases: 

-  $\xi > \langle\epsilon_f\rangle$

$$\begin{align}
I^{wt} &= \frac{\gamma}{\delta}\frac{\alpha\phi^{wt} - \delta}{\alpha\phi^{wt}}\frac{1-\langle\epsilon_b\rangle\xi}{1-\langle\epsilon_b\rangle \langle\epsilon_f\rangle}\\
I^{m} &= \frac{\gamma}{\delta}\frac{\alpha\phi^{m} - \delta}{\alpha\phi^{m}}\frac{1-\langle\epsilon_f\rangle\xi^{-1}}{1-\langle\epsilon_b\rangle\langle\epsilon_f\rangle}\\
\end{align}$$



"""

# ╔═╡ 4627e623-9e0b-4af2-92b3-84a0007f8be8
begin
	f_average = (M-1)/M*f + f_special/M
	b_average = (M-1)/M*b + b_special/M
end

# ╔═╡ 9b771a1e-f823-4b1c-a694-6e1d2d1fcb66
begin
	I_eq_wt_clarge = γ/δ * (α*ϕ_wt - δ)/(α*ϕ_wt) * 
		min(1, max(0, (1-b_average*ξ)/(1-b_average*f_average)))
	I_eq_m_clarge = γ/δ * (α*ϕ_m - δ)/(α*ϕ_m) * 
		min(1, max(0, (1-f_average/ξ)/(1-b_average*f_average)))

	freq_eq_clarge = I_eq_m_clarge / (I_eq_m_clarge + I_eq_wt_clarge)
end

# ╔═╡ 056fdf7a-8d16-4dd5-a8eb-eef280cbdc3b
md"### Final setup"

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
	v1 = PSS.Virus(; S=S0, I=I0, C=C0, ϕ=ϕ_wt) # wt virus
	v2 = PSS.Virus(; S=S0, I=0, C=0., ϕ=ϕ_m) # mutant: initially not here
	
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

# ╔═╡ af6249a4-1da6-4bff-8914-0b7ae9195e1b
ϕ_m_slider

# ╔═╡ 87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
p_freq = let
	# freq. plot
	tvals = range(sol.tspan..., length=100)
	f_R1 = PSS.frequency(sol, tvals, 1, 2)
	f_R2 = PSS.frequency(sol, tvals, 2, 2)
	f_R = PSS.frequency(sol, tvals, 2)

	lw = 2

	p = plot(
		legend = :topright,
		xlabel = "Time",
		ylabel = "",
		title = "Mutant frequency",
		ylim = (-0.03, 1.03),
	)
	
	plot!(tvals, f_R1, label="Region 1", linewidth=lw)
	plot!(tvals, f_R2, label="Other regions", linewidth=lw)
	plot!(tvals, f_R, label="Overall", linewidth=lw)
	hline!([freq_eq_c0]; label="c=0", line=(:black, lw+2, 0.3))
	hline!([freq_eq_clarge], label="c large", line=(:black, lw+1, .3, :dash))
	p
end

# ╔═╡ 95b3a0b8-d6fd-4554-89c3-73106c5db0d9
ϕ_m_slider

# ╔═╡ cc0ec860-c204-4b28-b4ce-d69da9116aee
Ms

# ╔═╡ 51dfe63a-b66e-4073-a531-dd64835158fe
Cs

# ╔═╡ 2e5f1f19-d5f2-4e54-8337-49db178b13c9
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

	hline!([I_eq_m_c0]; label="c=0", line=(:black, lw+2, 0.3))
	hline!([I_eq_m_clarge]; label="c large", line=(:dash, :black, lw+2, 0.3))

	hline!([I_eq_wt_c0]; label="c=0", line=(:blue, lw+2, 0.3))
	hline!([I_eq_wt_clarge]; label="c large", line=(:dash, :blue, lw+2, 0.3))
	p
end

# ╔═╡ add14358-6d60-47a5-876f-8e0651781c57
pS = let
	g = :S
	
	tvals = collect(range(sol.tspan..., length=400))
	lw = 2
	
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

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╟─db663438-8523-47ac-a783-88ff59a7ab1d
# ╟─c6c407d9-fea7-4e79-97ff-6acb88f4356d
# ╠═c79fe879-a32a-4391-8df1-5b9d9eeb75f5
# ╠═287520e8-97a3-485b-9634-fd2db7a2f4e6
# ╠═3074735c-5073-4a7d-8cc9-d983cc0268f2
# ╠═b3f03e99-9550-4b01-bda1-b7f9f6617c7d
# ╠═b14f025b-e95f-42fb-a321-4d914f6072ba
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═e38f872f-f9ec-45cf-953c-42a18aa0bab3
# ╠═644eed93-fdcc-479d-823d-076d40b517d5
# ╠═c2a159f4-cd25-4525-a270-2b3065a0dc0e
# ╟─bd7c5323-1763-411a-9d94-ba25d77ea8d9
# ╟─a7fc518a-f6d3-4e93-bac5-02d946c768f8
# ╠═c4da46f5-8178-4ae5-b865-db025a365da0
# ╟─81cbea6e-9a97-4e7e-9528-f0830c7ae58f
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╟─79bc278a-2e25-4ee4-b303-2f1f09cc3eac
# ╠═8ff3efd3-0771-4e84-b700-4177f90f72cb
# ╠═c26c16b7-af52-46ec-9013-b7c353876262
# ╠═f8025972-761d-4354-b688-35edaa639a80
# ╠═8e392fa9-4f6d-4b32-9c39-8e5c5000beb3
# ╠═abab9320-5e1f-475b-aa65-09c79f4cc284
# ╟─2bf1ad78-0389-445b-8334-c6f9bea86dcf
# ╟─8b265670-6152-409a-809a-f68c45203196
# ╠═9a14582a-1d87-49d6-8377-6df2a309da2d
# ╟─4643009b-aef8-469d-b25f-35df8b0593a9
# ╠═4627e623-9e0b-4af2-92b3-84a0007f8be8
# ╠═9b771a1e-f823-4b1c-a694-6e1d2d1fcb66
# ╟─056fdf7a-8d16-4dd5-a8eb-eef280cbdc3b
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
# ╠═af6249a4-1da6-4bff-8914-0b7ae9195e1b
# ╟─87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
# ╠═95b3a0b8-d6fd-4554-89c3-73106c5db0d9
# ╠═cc0ec860-c204-4b28-b4ce-d69da9116aee
# ╠═51dfe63a-b66e-4073-a531-dd64835158fe
# ╟─2e5f1f19-d5f2-4e54-8337-49db178b13c9
# ╟─add14358-6d60-47a5-876f-8e0651781c57
