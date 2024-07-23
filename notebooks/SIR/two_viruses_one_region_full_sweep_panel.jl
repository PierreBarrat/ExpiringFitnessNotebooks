### A Pluto.jl notebook ###
# v0.19.32

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

# ╔═╡ 2fe14f34-b81c-11ed-3729-67545b356a8b
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using Measures
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR # prefix functions with PSS, e.g. `PSS.Parameters`
end

# ╔═╡ b9fc4a25-66ae-4b70-86e4-633d62cd4b4e
include("$(homedir())/.julia/config/plot_defaults.jl")

# ╔═╡ 217462b9-a2b9-4f36-b1e6-37ea9bd58865
let
	def = pubfig(26)
	Plots.default(;def...)
end

# ╔═╡ 2754bfbd-9e4c-4881-a7a9-c5f4257141d1
PlutoUI.TableOfContents()

# ╔═╡ 1047ab91-66b9-4e95-8b1c-c98e5ee6d48d
md"# Parameters"

# ╔═╡ a451e995-e1e4-4845-a180-0072faf864ae
params = let
	N = 2
	M = 1
	α = 3
	γ = .002
	PSS.Parameters(; N, M, α, γ)
end

# ╔═╡ 8f30cd47-8dee-456c-ab67-366531762eff
# For direct access to α, ... 
@unpack α, γ, δ = params;

# ╔═╡ fb9110ed-b998-4d38-afe0-230c60b37872
md"""
Setting up cross-immunity: let us keep the equilibrium frequency $x$ constant and vary the distance $d$ of cross-immunity terms from $1$: 

$$\begin{aligned}
x &= \frac{1-f}{(1-f) + (1-b)} \;\; \text{(equilibrium frequency)}\\
& \\
d &= \frac{(1-b) + (1-f)}{2} \;\; \text{(distance from $1$ for $b$ and $f$)}\\ 
\end{aligned}$$

Inverting this gives 

$$f = 1-2xd, \; b = 1-2(1-x)d.$$
"""

# ╔═╡ 59a1579d-d0c1-4378-bb43-f85b568483d8
x_slider = @bind x Slider(0.05:0.05:0.95, show_value=true, default = 0.65)

# ╔═╡ d60e523f-be75-48d5-9f00-6ea0681f0a33
d_slider = @bind d Slider(0.025:0.025:0.95, show_value=true, default = 0.2)

# ╔═╡ 1d63d866-7ec0-4c59-bc50-a301a4324384
K = begin
	f = 1-2*x*d
	b = 1 - 2*(1-x)*d
	[1 b; f 1]
end

# ╔═╡ 199a949c-ef08-494c-b4d9-826129e1f92c
(f=f, b=b)

# ╔═╡ 0fb3d755-b409-47d0-bc4d-75a002ba96d6
md"""
**Note**: One can write `x_slider` or `d_slider` in any cell to tune the sliders
"""

# ╔═╡ b2148183-228f-45bf-a879-b5e6039ea5f4
md"# Simulation"

# ╔═╡ 7f8d658c-c329-4092-87ae-cfe33d44237e
md"Defining a region, initially with no mutant"

# ╔═╡ 55ffa5f0-0219-4b96-8b32-8251ab98b19e
I0 = 1e-6 # initial amount of mutant

# ╔═╡ 599065eb-0f6d-4d70-85ff-8e1ee6e2ad3f
region = let
	S0 = .4
	C0 = 1e-9
	R0 = 1 - S0 - I0 - C0
	wt = PSS.Virus(; S=S0, I=I0, C=C0, R=R0)
	mut = PSS.Virus(; S=1, I=0, C=0., R=0.)
	PSS.Region(; viruses=[wt, mut], K)
end

# ╔═╡ eb246d46-a899-4330-b430-a1c0ff5d8bc7
md"""
Simulating for time $T$ without the mutant to reach the initial equilibrium. Then, introduce a small quantity of the mutant in the region to obtain `state_init`. 
"""

# ╔═╡ 29993cd9-f004-49fa-88da-82b436373c6a
md"Finally, simulate again for time $T$."

# ╔═╡ 2e9996f8-093b-4d9d-b00b-b1d85c2a3267
md"Extract the final state, and also compute the analytical equilibrium."

# ╔═╡ b3a5ca74-0451-49dc-8057-e5994b87dce2
md"## Quick note on indexing"

# ╔═╡ 45894984-aff9-4a13-a81d-c6019bbd68fe
md"""
The type `PartialSweepSIR.SIRState`, like the variables `state_init` and `state_final`, can be indexed in the following way: 
- `state[i, a, X]` where `X` is in `[:S, :I, :C, :R]` returns the size of compartment `X` in region `i` for virus `a`. `i` and `a` can also be ranges, *e.g.* `state[1, 1:2, :I]`.
- `state[:, a, X]` returns a vector of size `M` with the size of compartment `X` in all regions $i\in[1,M]$
- `state[i, :, X]`: same, but for all viruses
"""

# ╔═╡ 22e56557-ec2c-45a3-b6cf-11aae0aa4975
md"""
The type `PartialSweepSIR.Solution`, like the variable `sol`, can be used like this: 
- `sol(t)` returns the state at time `t` in the form of an `SIRState` object. The paragraph above shows how to index it.
- `sol[tvals, i, a, X]` returns the vector of the size of compartment `X` in region `i` for virus `a` and for all values in `tvals`.  
"""

# ╔═╡ 91f7894d-5152-4239-8d61-77ed2fbcc818
md"# Figures"

# ╔═╡ d3584ba9-76a1-40e1-bdef-3121e2f5eba1
md"## Infectious"

# ╔═╡ bce748d4-7f16-47c3-9b22-cfdceb51f3a8
md"## Susceptibles"

# ╔═╡ 0a1f80c2-cf56-418e-b40d-d1645d11a8e9
md"## Frequency of mutant"

# ╔═╡ dd2038e4-0749-4c63-924c-540845439fe6
K

# ╔═╡ 99c10963-5219-4b49-9992-d9628dac69df
d_slider

# ╔═╡ 0ce1c6b9-b046-4488-b48d-af3c8422d6c7
x_slider

# ╔═╡ 0a922610-3be6-4c1f-b4f1-701cf6e76fcb
K

# ╔═╡ ae22a040-ce20-429b-9e5e-67208ab9a9db
init_growth_rate = let
	R0 = α/δ
	(1-f)*(R0-1) / (1 + f*(R0-1))
end

# ╔═╡ 203fb7e3-bfd2-4195-a6bc-f1f8452233f1
function logistic(t, x0, s)
	return exp(s*t) / (1/x0 - 1 + exp(s*t))
end

# ╔═╡ f98cf196-a1c1-4e02-a6ff-e1f2408084b5
let
	tvals = range(0, 200, 100)
	plot(tvals, logistic.(tvals, 1e-6, init_growth_rate))
end

# ╔═╡ 4363975c-c5e4-430e-8ea9-c9e2d51b36bf
logistic(0, 1e-8, init_growth_rate)

# ╔═╡ 8823aa48-bade-4f67-9b73-8816badff48c
blank_sp() = plot(legend=false,grid=false,foreground_color_subplot=:white)

# ╔═╡ 24cdb593-8465-4180-a42a-082bcb6a01a8
# Simulation time
T = 3/γ

# ╔═╡ 21040858-5c97-43a3-b38f-4e4ce5c90ef7
state_init = let
	state = PSS.SIRState(; regions=[region], parameters=params)
	sol = PSS.simulate(state, (0, T))
	PSS.set_infected(sol(T), 2, 1e-6) # infect with mutant
end;

# ╔═╡ 05f3e103-28a2-466e-829f-ea302f921712
state_eq = PSS.equilibrium(state_init); # analytical equilibrium

# ╔═╡ 8abbc0cd-e809-4d70-81d6-e4f783c94b53
sol = PSS.simulate(state_init, (0, T)); # simulated solution

# ╔═╡ 9d81f6f3-b341-4a56-ac8c-c63de91c541f
tvals = range(sol.tspan..., length=1_000) # relevant range of time values 

# ╔═╡ 003b59ba-412f-4dc7-8b77-5e8761b02ed1
pI = let
	v(a) = a == 1 ? "w.t." : "Mutant"
	
	# I 
	p = plot(
		yscale=:log10, 
		xlabel="time", 
		ylabel="",
		title="Infectious",
		legend=:bottomright,
		ylim = (1e-10, 1), 
		yticks = 10. .^ collect(-12:3:0)
	)

	for a in 1:2
		I = sol[tvals, 1, a, :I]
		plot!(tvals, I, color=a, label=v(a))
	end

	p
end

# ╔═╡ 6f0479e9-aae3-4052-911a-b5877d688c81
pS = let
	v(a) = a == 1 ? "w.t." : "Mutant"
	
	# I 
	p = plot(
		yscale=:linear, 
		xlabel="time", 
		ylabel="",
		title="Susceptibles",
		legend=:bottomright,
		xlim = (-10, 500)
	)

	for a in 1:2
		S = sol[tvals, 1, a, :S]
		plot!(tvals, S, color=a, label=v(a))
	end

	hline!([δ/α], line=(:grey, 10, .5), label="")
	
	p
end

# ╔═╡ 35d7f467-2a87-40c8-9142-993cf9057529
p_freq = let
	
	p = plot(
		# xlabel="time", 
		ylabel="",
		title="Mutant frequency",
		legend=:bottomright,
		frame = :box,
	)

	Iwt = sol[tvals, 1, 1, :I]
	Im = sol[tvals, 1, 2, :I]
	freq = Im ./ (Iwt + Im)

	plot!(tvals, freq, label="")
	hline!([x], line=(5, :grey, 0.5), label="")

	p
end

# ╔═╡ 5398dc6c-5847-43a7-b0ef-2d4f7353a1f3
state_final = sol(T);

# ╔═╡ c12b014b-5bce-4798-9f05-c486b0eb0715
p_freq_init = let
	tvals = range(0, T/15, length=1_000) # relevant range of time values 
	
	p = plot(
		xlabel="time", 
		ylabel="",
		# title="Mutant frequency",
		legend=:bottomright,
		frame = :box,
		# ylim = (1e-8, 1),
		# yscale=:log10,
		size = (900, 900)
	)

	Iwt = sol[tvals, 1, 1, :I]
	Im = sol[tvals, 1, 2, :I]
	freq = Im ./ (Iwt + Im)
	f0 = Im[1] / (Im[1] + Iwt[1])

	plot!(tvals, freq, label="")
	plot!(
		tvals, logistic.(tvals, f0, init_growth_rate);
		label="", line = (:black, :dash, 5, .5)
	)
	# hline!([x], line=(5, :grey, 0.5), label="")

	p
	
end

# ╔═╡ f3b71945-223e-4fc9-8cae-1a6c0335c02d
let
	panel = plot(
		pS, pI, p_freq;
		layout = grid(1,3), size = (3*900, 900), margin=10mm,
	)

	# savefig("$(homedir())/Documents/BaleLabo/Notes/ExpiringFitness/figures/two_strains_full_sweep_oscillations.png")

	panel
end

# ╔═╡ c86ef99b-ba3d-4540-837c-883c8943cef2
let
	panel = plot(
		pS, pI, p_freq, blank_sp(), blank_sp(), p_freq_init;
		layout = grid(2,3), size = (3*900, 2*900), margin=10mm,
	)

	savefig("$(homedir())/Documents/BaleLabo/Notes/ExpiringFitness/figures/two_strains_full_sweep_oscillations.png")

	panel
end

# ╔═╡ Cell order:
# ╠═2fe14f34-b81c-11ed-3729-67545b356a8b
# ╠═b9fc4a25-66ae-4b70-86e4-633d62cd4b4e
# ╠═217462b9-a2b9-4f36-b1e6-37ea9bd58865
# ╠═2754bfbd-9e4c-4881-a7a9-c5f4257141d1
# ╟─1047ab91-66b9-4e95-8b1c-c98e5ee6d48d
# ╠═a451e995-e1e4-4845-a180-0072faf864ae
# ╠═8f30cd47-8dee-456c-ab67-366531762eff
# ╟─fb9110ed-b998-4d38-afe0-230c60b37872
# ╠═59a1579d-d0c1-4378-bb43-f85b568483d8
# ╠═d60e523f-be75-48d5-9f00-6ea0681f0a33
# ╠═1d63d866-7ec0-4c59-bc50-a301a4324384
# ╠═199a949c-ef08-494c-b4d9-826129e1f92c
# ╟─0fb3d755-b409-47d0-bc4d-75a002ba96d6
# ╟─b2148183-228f-45bf-a879-b5e6039ea5f4
# ╟─7f8d658c-c329-4092-87ae-cfe33d44237e
# ╠═55ffa5f0-0219-4b96-8b32-8251ab98b19e
# ╠═599065eb-0f6d-4d70-85ff-8e1ee6e2ad3f
# ╟─eb246d46-a899-4330-b430-a1c0ff5d8bc7
# ╠═21040858-5c97-43a3-b38f-4e4ce5c90ef7
# ╟─29993cd9-f004-49fa-88da-82b436373c6a
# ╠═8abbc0cd-e809-4d70-81d6-e4f783c94b53
# ╟─2e9996f8-093b-4d9d-b00b-b1d85c2a3267
# ╠═5398dc6c-5847-43a7-b0ef-2d4f7353a1f3
# ╠═05f3e103-28a2-466e-829f-ea302f921712
# ╟─b3a5ca74-0451-49dc-8057-e5994b87dce2
# ╟─45894984-aff9-4a13-a81d-c6019bbd68fe
# ╟─22e56557-ec2c-45a3-b6cf-11aae0aa4975
# ╟─91f7894d-5152-4239-8d61-77ed2fbcc818
# ╠═9d81f6f3-b341-4a56-ac8c-c63de91c541f
# ╟─d3584ba9-76a1-40e1-bdef-3121e2f5eba1
# ╠═003b59ba-412f-4dc7-8b77-5e8761b02ed1
# ╟─bce748d4-7f16-47c3-9b22-cfdceb51f3a8
# ╠═6f0479e9-aae3-4052-911a-b5877d688c81
# ╟─0a1f80c2-cf56-418e-b40d-d1645d11a8e9
# ╠═dd2038e4-0749-4c63-924c-540845439fe6
# ╠═99c10963-5219-4b49-9992-d9628dac69df
# ╠═0ce1c6b9-b046-4488-b48d-af3c8422d6c7
# ╠═0a922610-3be6-4c1f-b4f1-701cf6e76fcb
# ╠═ae22a040-ce20-429b-9e5e-67208ab9a9db
# ╠═203fb7e3-bfd2-4195-a6bc-f1f8452233f1
# ╠═f98cf196-a1c1-4e02-a6ff-e1f2408084b5
# ╠═4363975c-c5e4-430e-8ea9-c9e2d51b36bf
# ╠═c12b014b-5bce-4798-9f05-c486b0eb0715
# ╠═8823aa48-bade-4f67-9b73-8816badff48c
# ╠═35d7f467-2a87-40c8-9142-993cf9057529
# ╠═24cdb593-8465-4180-a42a-082bcb6a01a8
# ╠═f3b71945-223e-4fc9-8cae-1a6c0335c02d
# ╠═c86ef99b-ba3d-4540-837c-883c8943cef2
