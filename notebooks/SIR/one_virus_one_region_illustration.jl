### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ c9b07ace-8157-11ee-2816-1f8500d98b5f
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR # prefix functions with PSS, e.g. `PSS.Parameters`
end

# ╔═╡ 82894038-f58a-4de3-af60-3bbae12a408b
params = let
	N = 1
	M = 1
	α = 3
	γ = .05
	PSS.Parameters(; N, M, α, γ)
end

# ╔═╡ b683f943-de1d-41fe-8a18-d91effd47545
# For direct access to α, ... 
@unpack α, γ, δ = params;

# ╔═╡ 043dbdeb-90cc-4f45-87ca-557ad69a90fb
T = 5/γ

# ╔═╡ 2d9ddc41-376a-4978-9729-312c6d8eb176
K = reshape([1.], 1, 1)

# ╔═╡ 6701cf01-665d-4a48-aaca-b2e071e2955c
region = let
	S0 = 1
	I0 = 1e-4
	C0 = 0
	R0 = 1 - S0 - I0 - C0
	wt = PSS.Virus(; S=S0, I=I0, C=C0, R=R0)
	PSS.Region(; viruses=[wt], K)
end

# ╔═╡ b7eea17a-91d1-4417-907d-f165ce7c39e6
state_init =  PSS.SIRState(; regions=[region], parameters=params)

# ╔═╡ ac5bfa98-bbfe-4773-ad36-7fbfffe7d6b1
sol = PSS.simulate(state_init, (0, T)); # simulated solution

# ╔═╡ c3d05a1f-c1e8-44d4-a37f-cbf898c77bb4
tvals = range(sol.tspan..., length=1_000) # relevant range of time values 

# ╔═╡ 0c8a733e-00f5-4ffc-a53c-cb4b723ae904
let
	lw = 3
	
	I = sol[tvals, 1, 1, :I]
	S = sol[tvals, 1, 1, :S]
	
	p = plot(
		xlabel="time", 
		ylabel="",
		title="",
		legend=:bottomright,
		frame=:box,
		legendfontsize=12,
		yscale=:log,
	)
	plot!(I, label="I"; lw)
	plot!(S, label="S"; lw)

	hline!([γ/δ*(1-δ/α)], label="", color=1, line=(:dash))

	# savefig("/home/pierrebc/Documents/BaleLabo/Slides/LJP_2023/SI_onestrain.png")
end

# ╔═╡ 12b3bdb1-21c6-4521-87fa-dc0acce2127f


# ╔═╡ Cell order:
# ╠═c9b07ace-8157-11ee-2816-1f8500d98b5f
# ╠═82894038-f58a-4de3-af60-3bbae12a408b
# ╠═b683f943-de1d-41fe-8a18-d91effd47545
# ╠═043dbdeb-90cc-4f45-87ca-557ad69a90fb
# ╠═2d9ddc41-376a-4978-9729-312c6d8eb176
# ╠═6701cf01-665d-4a48-aaca-b2e071e2955c
# ╠═b7eea17a-91d1-4417-907d-f165ce7c39e6
# ╠═ac5bfa98-bbfe-4773-ad36-7fbfffe7d6b1
# ╠═c3d05a1f-c1e8-44d4-a37f-cbf898c77bb4
# ╠═0c8a733e-00f5-4ffc-a53c-cb4b723ae904
# ╠═12b3bdb1-21c6-4521-87fa-dc0acce2127f
