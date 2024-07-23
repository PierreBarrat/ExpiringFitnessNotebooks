### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ bb54e23c-ee66-11ed-3ac8-63f12547877d
begin
	using Pkg; Pkg.activate("../../../")
	using Chain
	using CSV
	using DataFrames
	using JSON3
	using Measures
	using Plots
	using StatsBase
	using SpecialFunctions
end

# ╔═╡ cc2f13df-669f-4858-8096-ff6ddf0694d4
datdir = "data_trajectories_random_beta.jl/"

# ╔═╡ aa9a54e8-151f-453d-b313-35b1b5316af8
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 02186274-9f69-4371-8d9d-3815c0dcf890
Δt = 1

# ╔═╡ 77fcb250-3ea4-4435-865c-c30807105fee
dat = let
	df = DataFrame(values(files))
	df.diversity_strict = map(df.diversity_file) do f
		CSV.read(datdir * f, DataFrame).varpos_strict
	end
	df.diversity_005 = map(df.diversity_file) do f
		CSV.read(datdir * f, DataFrame).varpos_5
	end
	df[df.Δt .== Δt, :]
end;

# ╔═╡ db893776-c94c-41ad-ba94-a1b665a88622
αvals = dat.α |> unique |> sort

# ╔═╡ ac753dd2-91fa-43e8-8ea6-47c4c07c5aa8
mβvals = dat.mβ |> unique |> sort

# ╔═╡ 1393e551-5253-4ffc-9ad1-65d7e9695427
let
	r = first(eachrow(dat))
	names(r)
	# r.L * r.β2
	
end

# ╔═╡ 680a611a-89f0-4f91-964f-2752ab4e220f
md"""
Probability that $s = -\alpha\log(1-\beta) > \rho$, *i.e.* that sweeps do not overlap.

$$\begin{align*}
\mathbb{P} (\beta > 1-e^{-\rho/\alpha}) &\equiv \mathbb{P} (\beta > y)\\
&= \int_y^1 P(\beta)\text{d}\beta \\
&= 1 - \int_0^y B(x; a,b)\text{d}\beta\\
&= 1- I_y(a,b)
\end{align*}$$

where $B(x; a,b)$ is the Beta distribution and $I_y(a,b)$ the normalized incomplete Beta function. $a$ and $b$ are related to the first and second moments of $\beta$.
"""

# ╔═╡ 7e49eb53-37fd-409d-806a-0734762bbf80
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return a,b
end

# ╔═╡ 2c6f3b4f-ef9a-4214-9b17-33ce002007b8
function frac_overlap(mβ, β2, α, ρ)
	y = 1-exp(-ρ/α)
	a, b = get_β_distribution(mβ, β2)
	return beta_inc(a, b, y)[1]
end

# ╔═╡ 92c24bc7-375e-4049-bb23-73811a8dd7cd
function plot_div(row)
	p = plot(xlabel = "time", ylim = (-0.05, 1.05))
	tvals = range(0; length=length(row.diversity_strict), step=Δt)
	plot!(p, tvals, row.diversity_strict, label="Strict", linewidth=3)
	plot!(p, tvals, row.diversity_005, label="5%", linewidth=3)

	mβ = round(row.mβ, sigdigits=2)
	β2L = round(row.L * row.β2, sigdigits=2)
	p_noOVL = round(frac_overlap(row.mβ, row.β2, row.α, row.ρ); sigdigits=2)
	
	plot!(title = "<β> = $(mβ), <β2> = $(row.β2), α = $(row.α)\n <β2>L = $(β2L), Frac. Overlap = $(p_noOVL)")
	p
end

# ╔═╡ 4778f2e6-65b6-498c-b6d0-809de2790422
plt_dict = map(eachrow(dat)) do r
	(r.α, r.mβ) => plot_div(r)
end |> Dict

# ╔═╡ a8c8856f-1040-4fd6-9055-84b45d4e517f
p1 = let
	plot(
		[plt_dict[α, mβ] for α in αvals for mβ in mβvals]...,
		layout = grid(length(αvals), length(mβvals)),
		size = (length(mβvals) * 600, length(αvals) * 600),
		margin=10mm,
		titlefontsize=16,
		legendfontsize=16,
		tickfontsize=12
	)
end

# ╔═╡ a80d783b-2309-4a76-84a3-14903e3a41e0
savefig(p1, "random_beta_figures/diversity.png")

# ╔═╡ 3810d508-dcdd-4823-bb48-417f0a6bc2e6
md"# Tests"

# ╔═╡ 7c2dc726-89e2-4b45-8bfc-d91ab5b0c96c
let
	r = first(eachrow(dat))
	@info r.ρ
	ρ = r.ρ
	β2 = r.β2
	βs = range(mβvals[1], mβvals[end], 100)
	p = plot()
	for α in αvals
		prob = map(β -> frac_overlap(β, β2, α, ρ), βs)
		plot!(βs, prob, label="α=$α")
	end
	p
end

# ╔═╡ Cell order:
# ╠═bb54e23c-ee66-11ed-3ac8-63f12547877d
# ╠═cc2f13df-669f-4858-8096-ff6ddf0694d4
# ╠═aa9a54e8-151f-453d-b313-35b1b5316af8
# ╠═02186274-9f69-4371-8d9d-3815c0dcf890
# ╠═77fcb250-3ea4-4435-865c-c30807105fee
# ╠═db893776-c94c-41ad-ba94-a1b665a88622
# ╠═ac753dd2-91fa-43e8-8ea6-47c4c07c5aa8
# ╠═1393e551-5253-4ffc-9ad1-65d7e9695427
# ╠═4778f2e6-65b6-498c-b6d0-809de2790422
# ╠═a8c8856f-1040-4fd6-9055-84b45d4e517f
# ╠═a80d783b-2309-4a76-84a3-14903e3a41e0
# ╠═92c24bc7-375e-4049-bb23-73811a8dd7cd
# ╟─680a611a-89f0-4f91-964f-2752ab4e220f
# ╠═2c6f3b4f-ef9a-4214-9b17-33ce002007b8
# ╠═7e49eb53-37fd-409d-806a-0734762bbf80
# ╟─3810d508-dcdd-4823-bb48-417f0a6bc2e6
# ╠═7c2dc726-89e2-4b45-8bfc-d91ab5b0c96c
