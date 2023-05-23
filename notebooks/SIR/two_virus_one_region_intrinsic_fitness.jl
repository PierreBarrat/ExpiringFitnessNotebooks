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

# ╔═╡ 3f2181d2-b4ed-11ed-3e72-d7c1780d24e7
begin
	using Revise
	using Pkg; Pkg.activate("../../")
	using Parameters
	using PartialSweepSIR
	using PlutoUI
	using Plots
end

# ╔═╡ 1580d2cc-0654-4f93-affc-16da7f45220e
PlutoUI.TableOfContents()

# ╔═╡ 41076bec-3a0d-421e-a422-4b12ebe02b56
md"# Parameters"

# ╔═╡ c3a39edf-4727-485b-8382-aa21d071e3c7
params = let
	N = 2
	M = 1
	α = 3
	γ = .005
	PSS.Parameters(; N, M, α, γ)
end;

# ╔═╡ 88681913-11d0-47b1-b976-36219dcec855
@unpack α, γ, δ = params;

# ╔═╡ 5748f2c2-dd70-4aad-ba96-5813daef2f27
T = 10/γ

# ╔═╡ 35c68416-a19a-48f9-8957-395c31955230
md"## Cross-immunity"

# ╔═╡ 0b4f1c65-0ec3-45d9-a233-0bcaf8e1b552
md"""
Setting up cross-immunity; we control it using two parameters: 

$$\begin{aligned}
x &= \frac{b+f}{2} \;\; \text{(mean CI)}\\ 
&\\
Δ &= \frac{b - f}{2} \;\; \text{(asymetry of CI)}\\
\end{aligned}$$

Inverting this gives 

$$f = x - Δ, \; b = x + Δ.$$
"""

# ╔═╡ 1164c927-65ac-4a19-a12b-b392c81f497b
x_slider = let
	_xs = @bind x Slider(0.025:0.025:0.975, show_value=true, default = 0.8)
	md" $x = (b+f)/2$ = $(_xs)"
end

# ╔═╡ 57db3f6d-f88c-46dd-b966-80c40600617c
Δ_slider = begin
	dvals = vcat(
		range(-min(x, 1-x), 0; length=10), range(0, min(x, 1-x); length=10)
	) |> sort |> unique
	_ds = @bind Δ Slider(dvals, show_value=true, default = 0.)
	md" $\Delta = (b-f)/2$ = $(_ds)"
end

# ╔═╡ 19b876e1-f87f-4bed-96ca-4902f031a9cf
K = begin
	f = x - Δ
	b = x + Δ
	[1 b; f 1]
end

# ╔═╡ 6ed5662b-3286-4ca1-a559-082022b0d882
md"## Fitness"

# ╔═╡ 9fa67b59-6837-4dc5-879b-6c0060c79fd0
ϕ_wt = 1

# ╔═╡ 77a64be2-b40e-4f64-b57c-17e5a8cf6836
md"""
Co-existence condition: 

$$f < \xi \;\; \text{and} \;\;  b < \xi^{-1}$$

where

$$\xi = \frac{\alpha\phi^m - \delta}{\alpha\phi_{wt} - \delta}.$$
"""

# ╔═╡ e7dfd56e-9719-4f9b-a988-e753e118d54f
ϕs = let
	_slider = @bind ϕ_m Slider(0.5:0.01:2, show_value=true, default=1)
	_ξ = (ϕ_m - δ/α)/(ϕ_wt - δ/α)
	ϕs = md"""
	 $\phi_m$  = $(_slider) 
	"""
end

# ╔═╡ 37a01f00-b8b6-49cc-8ba3-c36593f9fd7d
ξs = md" $\xi$ = $((ϕ_m - δ/α)/(ϕ_wt - δ/α))"

# ╔═╡ 6dff5b25-43a0-46b7-84bd-8cf438542511
ϕ_m_slider = (ϕs, ξs)

# ╔═╡ 69171b6c-0d3c-43ed-b52e-1e4da14ae33b
md"# Simulation"

# ╔═╡ 3d1dcb13-ac89-46fd-9466-2f013aaa2935
region = let
	S0 = .4
	I0 = 1e-6
	C0 = 1e-9
	R0 = 1 - S0 - I0 - C0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0, R=R0, ϕ = ϕ_wt)
	v2 = PSS.Virus(; S=1, I=0, C=0., R=0., ϕ = ϕ_m)
	PSS.Region(; viruses=[v1,v2], K)
end

# ╔═╡ 9bc81dc4-7a27-431f-876c-afaf11b3f54e
state_init = let
	state = PSS.SIRState(; regions=[region], parameters=params)
	sol = PSS.simulate(state, (0, T))
	state_init = PSS.set_infected(sol(T), 2, 1e-6) # infect with mutant
	state_init
end;

# ╔═╡ 2cd088e3-644a-45bd-b9e2-228063501f59
sol = PSS.simulate(state_init, (0, T)); # simulated solution

# ╔═╡ 2e822177-bab4-4242-8d55-cdb43e1da60e
state_final = sol(T);

# ╔═╡ e3d1ce2b-85a8-4f5f-b58f-acc142b18370
md"# Analytical equilibrium"

# ╔═╡ a26a326d-112d-4483-a27a-e2974971b017
begin
	Rwt = α/δ*ϕ_wt
	Rm = α/δ*ϕ_m
	ξ = (ϕ_m - δ/α)/(ϕ_wt - δ/α)

	I_eq_wt = γ/δ * (Rwt-1)/Rwt * min(1, max(0, (1 - b*ξ)/(1-b*f)))
	I_eq_m = γ/δ * (Rm-1)/Rm * min(1, max(0, (1-f/ξ)/(1-b*f)))
end

# ╔═╡ 21a380ab-de02-4407-8d6c-cafe8a85f5d1
md"Checking that numerically inverting the matrix and the analytical formula give the same result."

# ╔═╡ b99267e0-f2ba-4c5c-bae1-b1107b28422a
state_eq = PSS.equilibrium(state_init); # inverts K numerically

# ╔═╡ 5d6e82ce-f814-4527-8d88-5c3d53cb0b0d
let
	@assert isapprox(I_eq_m, state_eq[1, 2, :I])
	@assert isapprox(I_eq_wt, state_eq[1, 1, :I])
end

# ╔═╡ a0ed2aa2-a275-4aa6-b3fe-2b0260f9968e
state_eq[1, 2, :I]

# ╔═╡ 123ac0a0-c55d-44e2-bb8b-8b1643154c21
I_eq_m

# ╔═╡ 991ea65d-346d-4704-8315-56ed45547e6a
f_eq = (ξ - f) / (ξ - f + ϕ_m/ϕ_wt - b*ξ*ϕ_m/ϕ_wt)

# ╔═╡ 6c558d3a-77be-4124-aa20-767c58510c91
md"# Figures"

# ╔═╡ d427794f-675c-4a16-8b31-5f347fca92f8
ϕ_m_slider

# ╔═╡ c4f5530c-4646-4695-9b32-41fc9a36190b
pI = let
	tvals = range(sol.tspan..., length=200)
	v(a) = a == 1 ? "w.t." : "Mutant"
	# # I 
	p = plot(
		yscale=:log10, 
		ylim = (1e-7, 5e-2), 
		xlabel="", 
		ylabel="",
		# title="Infectious",
		legend=:bottomright,
	)
	g = :I
	for a in 1:2
		X = sol[tvals, 1, a, g]
		plot!(tvals, X, color=a, label=v(a))
	end

	hline!([I_eq_m], label="", color=2, line=(:dash), alpha=0.5)
	hline!([I_eq_wt], label="", color=1, line=(:dash), alpha=0.5)

	# p
end

# ╔═╡ f092b2d0-4914-4dc6-ab6a-1b1337d8f36e
ϕ_m_slider

# ╔═╡ 193be30d-e9c9-475d-9b7e-8e9007921e6a
Δ_slider

# ╔═╡ 743382c9-5bf7-4cef-8afc-7b47b89385f2
x_slider

# ╔═╡ e2bcf2fb-285a-47c7-8464-9855608af011
let
	tvals = range(sol.tspan..., length=200)
	lw = 2
	v(a) = a == 1 ? "w.t." : "Mutant"
	
	# I 
	p = plot(
		ylim = (-0.05, 1.05), 
		xlabel="", 
		ylabel="",
		title="Frequency",
		legend=:bottomright,
	)
	g = :I
	X = sol[tvals, 1, 2, g] ./ (sol[tvals, 1, 1, g] + sol[tvals, 1, 2, g])
	plot!(tvals, X, color=1, label="", line=(lw))

	hline!([f_eq], label="", line=(:dash, :black, lw, 0.5))

	p
end

# ╔═╡ Cell order:
# ╠═3f2181d2-b4ed-11ed-3e72-d7c1780d24e7
# ╠═1580d2cc-0654-4f93-affc-16da7f45220e
# ╟─41076bec-3a0d-421e-a422-4b12ebe02b56
# ╠═c3a39edf-4727-485b-8382-aa21d071e3c7
# ╠═88681913-11d0-47b1-b976-36219dcec855
# ╠═5748f2c2-dd70-4aad-ba96-5813daef2f27
# ╟─35c68416-a19a-48f9-8957-395c31955230
# ╟─0b4f1c65-0ec3-45d9-a233-0bcaf8e1b552
# ╠═1164c927-65ac-4a19-a12b-b392c81f497b
# ╠═57db3f6d-f88c-46dd-b966-80c40600617c
# ╠═19b876e1-f87f-4bed-96ca-4902f031a9cf
# ╟─6ed5662b-3286-4ca1-a559-082022b0d882
# ╠═9fa67b59-6837-4dc5-879b-6c0060c79fd0
# ╟─77a64be2-b40e-4f64-b57c-17e5a8cf6836
# ╠═e7dfd56e-9719-4f9b-a988-e753e118d54f
# ╠═37a01f00-b8b6-49cc-8ba3-c36593f9fd7d
# ╠═6dff5b25-43a0-46b7-84bd-8cf438542511
# ╟─69171b6c-0d3c-43ed-b52e-1e4da14ae33b
# ╠═3d1dcb13-ac89-46fd-9466-2f013aaa2935
# ╠═9bc81dc4-7a27-431f-876c-afaf11b3f54e
# ╠═2cd088e3-644a-45bd-b9e2-228063501f59
# ╠═2e822177-bab4-4242-8d55-cdb43e1da60e
# ╟─e3d1ce2b-85a8-4f5f-b58f-acc142b18370
# ╠═a26a326d-112d-4483-a27a-e2974971b017
# ╟─21a380ab-de02-4407-8d6c-cafe8a85f5d1
# ╠═b99267e0-f2ba-4c5c-bae1-b1107b28422a
# ╠═5d6e82ce-f814-4527-8d88-5c3d53cb0b0d
# ╠═a0ed2aa2-a275-4aa6-b3fe-2b0260f9968e
# ╠═123ac0a0-c55d-44e2-bb8b-8b1643154c21
# ╠═991ea65d-346d-4704-8315-56ed45547e6a
# ╟─6c558d3a-77be-4124-aa20-767c58510c91
# ╠═d427794f-675c-4a16-8b31-5f347fca92f8
# ╠═c4f5530c-4646-4695-9b32-41fc9a36190b
# ╠═f092b2d0-4914-4dc6-ab6a-1b1337d8f36e
# ╠═193be30d-e9c9-475d-9b7e-8e9007921e6a
# ╠═743382c9-5bf7-4cef-8afc-7b47b89385f2
# ╟─e2bcf2fb-285a-47c7-8464-9855608af011
