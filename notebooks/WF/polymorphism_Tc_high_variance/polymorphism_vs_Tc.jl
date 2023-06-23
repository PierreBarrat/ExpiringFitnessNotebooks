### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 5127cf4e-eb40-11ed-0a15-53a1acf43aef
begin
	using Revise
	using Pkg; Pkg.activate("../../../")
	using Chain
	using CSV
	using DataFrames
	using DelimitedFiles
	using Distributions
	using FrequencyTrajectories
	using JSON3
	using Measures
	using Plots
	using Random
	using StatsBase
end

# ╔═╡ 5ffd9c4e-4807-43aa-943a-4c8c3f9b022a
datdir = "data_trajectories_random_beta.jl/"

# ╔═╡ 2713455c-df95-49c7-925a-9d92b91719a5
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 9ac16620-e45c-4218-9d7f-c7cca89e66e1
md"# Functions"

# ╔═╡ 9ff2d96e-4637-4699-a4a2-47f4d9887fa5
"""
	get_β_distribution(mβ, β2)

Return a Beta distribution with mean `mβ` and second moment `β2`.
"""
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return Beta(a,b)
end

# ╔═╡ ed80a931-f103-4a28-9e17-17ec8ed5aed9
function sfs(freq, positions)
	return map(eachrow(freq)) do f
		2 * mean(x->x*(1-x), f[positions[1:2:end]])
	end
end

# ╔═╡ c9a65c7c-2b3b-403a-9141-6bc9fb316512
function fraction_covered_by_sweeps(mβ, β2, ρ, α)
	B = get_β_distribution(mβ, β2)
	return - 3 * ρ / α / mean(β -> log(1-β), rand(B, 1000))
end

# ╔═╡ 10987b07-70e0-4a19-9173-83cf487c08ff
dat = let
	df = DataFrame(values(files))
	df.poly = map(df.diversity_file) do f
		X = CSV.read(datdir * f, DataFrame)
		X.polymorphism / 2
	end
	df.tvals =  map(df.diversity_file) do f
		X = CSV.read(datdir * f, DataFrame)
		X.t
	end

	df.sweep_fraction = map(eachrow(df)) do r
		fraction_covered_by_sweeps(r.mβ, r.β2, r.ρ, r.α)
	end
	df.Tc_theory = map(r -> 1/r.ρ/r.β2, eachrow(df))
	df.Tc_measured = map(eachrow(df)) do r
		x = mean(r.poly[1:end])
		x/r.μ_neutral
	end

	# df = df[df.ρ .== minimum(df.ρ) .&& df.mβ .== maximum(df.mβ), :]
	
	sort(df, [:ρ, :β2])
end;

# ╔═╡ 94813cd3-210a-4b77-aac7-f68e1dc3a1ba
begin
	# Constants
	α = dat.α[1]
	μ_neutral = dat.μ_neutral[1]
	L_sel = dat.L[1]
	L_neutral = dat.L_neutral[1]
	L = L_sel + L_neutral
	neutral_sites = (1+2*L_sel):2*L
end

# ╔═╡ 60ae9a79-0ec1-483e-ab4a-0477cb2d93cc
function estimate_Tc(mβ, β2)
	T1 = 1/ρ/β2

	βsamples = rand(get_β_distribution(mβ, β2), 1000)
	T2 = mean(β -> -3/α/log(1-β), βsamples)
	
	return T1 + T2, T1, T2
end

# ╔═╡ 3eb3f135-af7a-44d5-acce-64c371e13e84
dat[4,:].mβ^2

# ╔═╡ b2a91a71-d03f-41e4-8d1a-d58ce7a6a13f
dat[4,:].β2

# ╔═╡ 9b6ba4bb-a5e6-4f19-9da3-bb9690c6c0aa
parameters = map(r -> (mβ=r.mβ, β2=r.β2, ρ=r.ρ), eachrow(dat)) |> sort |> unique

# ╔═╡ 6a5d6c8f-618d-4cb4-a519-31042c459f5f
begin
	β2_vals = dat.β2 |> unique |> sort
	ρ_vals = dat.ρ |> unique |> sort
end

# ╔═╡ 982ff3e1-47f4-48de-b4fc-c696b239e6e1
let
	grouped_df = groupby(dat, :mβ, sort=true)
	p = plot()
	for (i, (k, df)) in enumerate(zip(keys(grouped_df), grouped_df))
		mβ = round(k.mβ, sigdigits=2)
		err = abs.((df.Tc_measured - df.Tc_theory) ./ df.Tc_theory)
		# err = (df.Tc_measured - df.Tc_theory) ./ df.Tc_measured
		scatter!(
			1 ./df.ρ, err, 
			label = "<β> = $(mβ)", color=i
		)
		# plot!(1 ./ df.ρ, 1 ./ df.ρ ./ df.β2, color=i, label="")
	end
	# plot!(scale=:log10)
	p
end

# ╔═╡ c938b1e5-d778-4b41-b789-1f03acac604a
let
	grouped_df = groupby(dat, :mβ, sort=true)
	p = plot()
	for (i, (k, df)) in enumerate(zip(keys(grouped_df), grouped_df))
		mβ = round(k.mβ, sigdigits=2)
		scatter!(1 ./df.ρ, df.Tc_measured, label = "<β> = $(mβ)", color=i)
		plot!(1 ./ df.ρ, 1 ./ df.ρ ./ df.β2, color=i, label="")
	end
	plot!(scale=:log10)
end

# ╔═╡ 83eb61b1-622b-4b1d-a1bf-3a7a0a750cb5
function sweep_time(mβ, β2, α)
	B = get_β_distribution(mβ, β2)
	return @chain begin
		rand(B, 10_000) # β sample
		map(β -> log(1-β), _)
		filter(isfinite, _)
		mean
		- 3 / α / _
	end
end

# ╔═╡ 04a1b0ee-233f-4c94-a150-c95ffd10c868
function plot_polymorphism(r)
	idx = 1:10:length(r.poly)
	plot(r.tvals[idx], r.poly[idx] / r.μ_neutral, label="", alpha=0.3)
	hline!([r.Tc_theory], label="Ne", line=(:red, 3, 0.9))

	hline!([r.Tc_measured], label="", color=1, line=(:dashdot, 4))
	
	st = round(sweep_time(r.mβ, r.β2, r.α), sigdigits=2)
	dt = round(Int, 1/r.ρ)
	plot!(
		xlabel="time",
		ylabel = "x(1-x)/μ",
		# title = "1/ρ = $(dt), sweep time ~ $(st)",
		title = "ρ=$(round(r.ρ, sigdigits=2)), β = $(r.mβ)\n 1/ρ = $(dt), sweep time ~ $(st)",
		# ylim = (-100, 6*Ne),
		tickfontsize = 6,
	)
end

# ╔═╡ bd8b7c82-c16f-4fbf-ac92-6e213580a641
plot_polymorphism(dat[4,:])

# ╔═╡ 21c7431a-3b6d-41f0-b5fb-e989bad2259a
plts = map(Iterators.product(ρ_vals, β2_vals)) do (ρ, β2)
	df = dat[dat.β2 .== β2 .&& dat.ρ .== ρ, :]
	plot_polymorphism(first(eachrow(df)))
end;

# ╔═╡ 88229ffb-bdb1-4691-9f1a-2e8eccec2954
let
	p = plot(
		plts...;
		layout = (length(β2_vals), length(ρ_vals)),
		size = (400*length(ρ_vals), 400*length(β2_vals)),
		margin = 5mm, 
		guidefontsize = 16,
		legendfontsize = 16,
		tickfontsize = 10, 
	)
	savefig("panel_diversity_v_time.png")
	p
end

# ╔═╡ 134f46e3-6c49-477c-909c-7dad0b450b0c
sweep_time(dat[3,:].mβ, dat[3,:].β2, dat[3,:].α)

# ╔═╡ Cell order:
# ╠═5127cf4e-eb40-11ed-0a15-53a1acf43aef
# ╠═5ffd9c4e-4807-43aa-943a-4c8c3f9b022a
# ╠═2713455c-df95-49c7-925a-9d92b91719a5
# ╠═10987b07-70e0-4a19-9173-83cf487c08ff
# ╠═94813cd3-210a-4b77-aac7-f68e1dc3a1ba
# ╠═04a1b0ee-233f-4c94-a150-c95ffd10c868
# ╠═bd8b7c82-c16f-4fbf-ac92-6e213580a641
# ╠═3eb3f135-af7a-44d5-acce-64c371e13e84
# ╠═b2a91a71-d03f-41e4-8d1a-d58ce7a6a13f
# ╠═9b6ba4bb-a5e6-4f19-9da3-bb9690c6c0aa
# ╠═6a5d6c8f-618d-4cb4-a519-31042c459f5f
# ╠═21c7431a-3b6d-41f0-b5fb-e989bad2259a
# ╠═88229ffb-bdb1-4691-9f1a-2e8eccec2954
# ╠═982ff3e1-47f4-48de-b4fc-c696b239e6e1
# ╠═c938b1e5-d778-4b41-b789-1f03acac604a
# ╠═134f46e3-6c49-477c-909c-7dad0b450b0c
# ╠═9ac16620-e45c-4218-9d7f-c7cca89e66e1
# ╠═60ae9a79-0ec1-483e-ab4a-0477cb2d93cc
# ╠═9ff2d96e-4637-4699-a4a2-47f4d9887fa5
# ╠═ed80a931-f103-4a28-9e17-17ec8ed5aed9
# ╠═c9a65c7c-2b3b-403a-9141-6bc9fb316512
# ╠═83eb61b1-622b-4b1d-a1bf-3a7a0a750cb5
