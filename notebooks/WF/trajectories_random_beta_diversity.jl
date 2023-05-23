### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ bb54e23c-ee66-11ed-3ac8-63f12547877d
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using JSON3
	using Measures
	using Plots
	using StatsBase
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

# ╔═╡ 92c24bc7-375e-4049-bb23-73811a8dd7cd
function plot_div(row)
	p = plot(xlabel = "time", ylim = (-0.05, 1.05))
	tvals = range(0; length=length(row.diversity_strict), step=Δt)
	plot!(p, tvals, row.diversity_strict, label="Strict", linewidth=3)
	plot!(p, tvals, row.diversity_005, label="5%", linewidth=3)

	plot!(title = "<β> = $(round(row.mβ, sigdigits=2)), <β2> = $(row.β2), α = $(row.α)")
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
		size = (length(mβvals) * 450, length(αvals) * 450),
		margin=10mm,
		titlefontsize=16,
		legendfontsize=16,
		tickfontsize=12
	)
end

# ╔═╡ a80d783b-2309-4a76-84a3-14903e3a41e0
savefig(p1, "random_beta_figures/diversity.png")

# ╔═╡ Cell order:
# ╠═bb54e23c-ee66-11ed-3ac8-63f12547877d
# ╠═cc2f13df-669f-4858-8096-ff6ddf0694d4
# ╠═aa9a54e8-151f-453d-b313-35b1b5316af8
# ╠═02186274-9f69-4371-8d9d-3815c0dcf890
# ╠═77fcb250-3ea4-4435-865c-c30807105fee
# ╠═db893776-c94c-41ad-ba94-a1b665a88622
# ╠═ac753dd2-91fa-43e8-8ea6-47c4c07c5aa8
# ╠═4778f2e6-65b6-498c-b6d0-809de2790422
# ╠═a8c8856f-1040-4fd6-9055-84b45d4e517f
# ╠═a80d783b-2309-4a76-84a3-14903e3a41e0
# ╠═92c24bc7-375e-4049-bb23-73811a8dd7cd
