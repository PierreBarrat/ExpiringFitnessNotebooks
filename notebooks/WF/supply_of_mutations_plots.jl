### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 99ff2586-cfb5-11ed-0639-7348bc0a211a
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

# ╔═╡ ae05c4ee-3f62-4728-a13c-8989beac81e4
datdir = "data_supply_of_mutations/"

# ╔═╡ d95b12c0-bfc6-46a0-9e40-07b6ca6a9874
filenames = JSON3.read(read(datdir * "files.json", String), Dict);

# ╔═╡ e4cc2b79-73cc-471d-90a4-6f56dde7cbd8
data = map(zip(keys(filenames), values(filenames))) do (f, x)
	dat = CSV.read(datdir * f, DataFrame)
	ρ = x["ρ"]
	β = x["β"]
	L = x["L"]
	α = x["α"]
	(ρ=ρ, β=β, L=L, α=α) => dat
end |> Dict;

# ╔═╡ 33b8e7cf-eccb-42a7-965e-dae5e3936f7c
ρ, β, _ = first(keys(data))

# ╔═╡ fe321e1a-5770-4c94-8446-f23103cd5cf0
dat = first(values(data));

# ╔═╡ b0d9b620-ae03-487f-b6f5-486f0961d0c1
ρvals = map(x -> x.ρ, collect(keys(data))) |> sort |> unique

# ╔═╡ af718c5c-b8dd-4d7b-80f3-66114ceebc9e
βvals = map(x -> x.β, collect(keys(data))) |> sort |> unique

# ╔═╡ 26b5bcef-f313-4991-b0b5-e5fb86c2fca2
L = map(x -> x.L, collect(keys(data))) |> sort |> unique |> first

# ╔═╡ d04bc7fc-fe79-4d03-af9d-fd5cb4499a9c
α = map(x -> x.α, collect(keys(data))) |> sort |> unique |> first

# ╔═╡ 060d74fb-ad57-4cd9-b08e-2f678a6e1f88
size = (600, 600)

# ╔═╡ 52c63e55-e412-4c62-9a7d-e8e067fa61d3
compute_s(β, α) = -α * log(1-β)

# ╔═╡ d26365ae-2c90-41b1-854f-cae4d60f9213
function do_plot(ρ, β, L, α, data)
	s = compute_s(β, α)
	v1 = round(L*β^2, sigdigits=2)
	v2 = round(s/ρ, sigdigits=2)
	
	p = plot()

	plot!(data.t, data.varpos_5, label = "thr: 0.05")
	plot!(data.t, data.varpos_strict, label = "thr: 1/N")
	hline!([β^-2/L], line=(:black, :dashdot, 3, 0.5), label = "1/β^2/L")
	vline!([1/ρ/β^2], line = (:black, :dash), label = "Ne")
	
	plot!(
		xlabel = "time",
		ylabel = "Nb var. positions (scaled)",
		title = "ρ=$ρ, β=$β, Lβ^2=$(v1), s/ρ=$(v2)",
		frame=:box,
		ylim = (-0.02, 1.02),
	)

	return p
end

# ╔═╡ 59d0cbf7-982a-4068-9ada-951350cfa61c
do_plot(ρ, β, L, α, dat)

# ╔═╡ 45d0f87d-ac7c-43df-894b-57425cfab058
let
	plts = [
		do_plot(ρ, β, L, α, data[(ρ=ρ, β=β, L=L, α=α)]) 
		for β in βvals, ρ in ρvals
	]
	p = plot(
		plts..., 
		layout = (length(ρvals), length(βvals)), 
		size = (3000, 2400),
		margin=8mm
	)

	savefig(p, datdir * "variable_positions_v_time.png")
	p
end

# ╔═╡ Cell order:
# ╠═99ff2586-cfb5-11ed-0639-7348bc0a211a
# ╠═ae05c4ee-3f62-4728-a13c-8989beac81e4
# ╠═d95b12c0-bfc6-46a0-9e40-07b6ca6a9874
# ╠═e4cc2b79-73cc-471d-90a4-6f56dde7cbd8
# ╠═33b8e7cf-eccb-42a7-965e-dae5e3936f7c
# ╠═fe321e1a-5770-4c94-8446-f23103cd5cf0
# ╠═59d0cbf7-982a-4068-9ada-951350cfa61c
# ╠═b0d9b620-ae03-487f-b6f5-486f0961d0c1
# ╟─af718c5c-b8dd-4d7b-80f3-66114ceebc9e
# ╠═26b5bcef-f313-4991-b0b5-e5fb86c2fca2
# ╠═d04bc7fc-fe79-4d03-af9d-fd5cb4499a9c
# ╠═45d0f87d-ac7c-43df-894b-57425cfab058
# ╠═060d74fb-ad57-4cd9-b08e-2f678a6e1f88
# ╠═d26365ae-2c90-41b1-854f-cae4d60f9213
# ╠═52c63e55-e412-4c62-9a7d-e8e067fa61d3
