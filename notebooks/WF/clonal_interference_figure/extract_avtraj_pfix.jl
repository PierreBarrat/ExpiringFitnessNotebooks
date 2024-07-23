### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 0aaf2e58-04af-4387-90a9-f1aa2632cd97
begin
	using Revise
	using Pkg; Pkg.activate("../../../")
	using Chain
	using CSV
	using DataFrames
	using DataFramesMeta
	using FrequencyTrajectories
	using JSON3
	using Measures
	using Plots
	using Random
	using StatsBase
end

# ╔═╡ cf871f3d-cbb2-49c6-ae88-81c0f3ec91ca
datdir = "../data_clonal_interference"

# ╔═╡ 6704b259-a757-4549-a1cd-f8a9c261b87a
files = open(joinpath(datdir, "files.json"), "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ bfb9d9d3-4560-4f49-86a0-01aafeaa37c8


# ╔═╡ 1eff07a7-c0a5-46f6-af57-9af2d45ddbd2


# ╔═╡ 52c4f0a4-c367-49b7-be89-5d038b5586aa
md"# Average trajectory"

# ╔═╡ 6d73d6bc-cb3d-488f-a47a-18aceba60dbf
md"# Probability of fixation"

# ╔═╡ 47164d93-e263-405c-a071-e32270b83fd8
md"## Functions"

# ╔═╡ 29c34604-4b4b-4cf7-be40-ea8054ca173b
function read_trajectories(file)
	return @chain file begin
		joinpath(datdir, _)
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories_from_dataframe
	end
end

# ╔═╡ cf9df8e6-cf64-412c-a263-ba746f3ff706
dat_all = let
	df = DataFrame(values(files))
	df.trajectories = map(f -> read_trajectories(f), df.trajectory_file)
	df
end;

# ╔═╡ 5b9c1d07-d54b-4c3e-8962-f97793595d82
begin
	svals = dat_all.s |> unique |> sort
	ρvals = dat_all.ρ |> unique |> sort
	
	fb = FrequencyBin(0.3, 0.05)

	cfs = "random"
end

# ╔═╡ 3ce576c6-ab52-40ff-abad-b59521fd75a3
dat = @chain dat_all begin
	@subset(:cfs .== cfs)
	sort!([:ρ])
end;

# ╔═╡ 1c390019-a2d4-4a2d-8db6-f0159a5ff82a
names(dat)

# ╔═╡ 3a0443be-f4f9-475c-b17a-2fe2ed32e7c4
always_below=true

# ╔═╡ 4daade11-9081-4c8c-a42a-7cda945c2611
function get_mean_trajectory(trajectories, fb; kwargs...)
	fmean, tmean, Zs = mean(trajectories, fb; always_below, kwargs...)
	return fmean, tmean
end

# ╔═╡ fba6648b-b279-48cd-a956-0631eafc6899
let 
	av_traj_df = select(
		dat, :ρ, :s, 
		[:trajectories] => ByRow(t -> get_mean_trajectory(t, fb; K=1/2)) => [:mf, :mt]
	)
	CSV.write("data_avtraj_$(cfs).csv", av_traj_df)
end

# ╔═╡ e29bf4c8-c686-4565-9c2a-33a1eb2b70c8
function _plot_trajectories(trajectories, fb)
	tmin = -400
	tmax = 1200
	# fb = FT.FrequencyBin(.4, 0.03)
	
	local p = plot(
		ylim=(-0.02, 1.02), 
		xlim = (tmin, tmax),
		xlabel = "Time",
		ylabel = "Frequency",
		size = (900, 900),
		frame = :box
	)
	
	local fmean, tmean, Zs = mean(trajectories, fb; always_below)
	for t in shuffle(FrequencyTrajectories.filter(trajectories, fb; always_below))
		plot!(t, fb, label="", alpha = 0.25,)
	end
	
	plot!(tmean, fmean, line=(:black, 3), label="")
	plot!(collect(xlims(p)), [0.95, 0.95], label="", line=(:black, :dash))
	plot!(collect(xlims(p)), [0.05, 0.05], label="", line=(:black, :dash))
	plot!([tmin, 0], fb)
	
	return p
end

# ╔═╡ c5a36cc9-5622-4abe-9dc4-5adbf0b999e9
function plot_trajectories(trajectories, params, fb)
	p = _plot_trajectories(trajectories, fb)
	# vline!([params.Ne], line=(:black, :dash), label="Ne=$(params.Ne)")
	plot!(
		p,
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), s/ρ = $(round(params.s/params.ρ; sigdigits=2)), Δt = $(params.Δt)"
	)
end

# ╔═╡ 6e9b4d7b-ff90-479b-aecf-a269f7716699
function pfix_v_f(trajectories)
	fbs = [FrequencyBin(f, 0.05) for f in .1:.1:.9]
	pfix = 	map(fbs) do fb
		fixation_probability(T -> FT.inbin(T, fb; always_below), trajectories)
	end
	return (fbs, pfix)
end

# ╔═╡ 197fe84b-7314-4baf-ae79-4985b218db62
let
	pfix_df = select(
		dat, :ρ, :s, [:trajectories] => ByRow(pfix_v_f) => [:fb, :pfix]
	) 
	select!(pfix_df, Not(:fb), :fb => ByRow(col -> map(fb -> fb.f, col)) => :f)
	CSV.write("data_pfix_$(cfs).csv", pfix_df)
end

# ╔═╡ fef63f90-b890-4a39-b577-aa420481b185
function plot_pfix(trajectories, params)
	fbs, pfix = pfix_v_f(trajectories)

	p = plot([fb.f for fb in fbs], pfix, label="", marker=:o)
	plot!([0,1], [0,1], line=(:black, :dash), label="")
	plot!(
		p,
		title = "s = $(params.s), ρ = $(round(params.ρ, sigdigits=2)), s/ρ = $(round(params.s/params.ρ; sigdigits=2)), Δt = $(params.Δt)"
	)
end

# ╔═╡ 96b32c90-2e8f-43d0-9407-3206f9ad9f0e
plt_dict_pfix = map(eachrow(dat)) do r
	(r.s, r.ρ) => plot_pfix(r.trajectories, r)
end |> Dict

# ╔═╡ a97405ed-6a50-4b97-9af2-5b8a074eabba
p2 = let
	plot(
		[plt_dict_pfix[s, ρ] for s in svals for ρ in ρvals]...,
		layout = grid(length(svals), length(ρvals)),
		size = (length(ρvals) * 600, length(svals) * 600)
	)
end

# ╔═╡ Cell order:
# ╠═0aaf2e58-04af-4387-90a9-f1aa2632cd97
# ╠═cf871f3d-cbb2-49c6-ae88-81c0f3ec91ca
# ╠═6704b259-a757-4549-a1cd-f8a9c261b87a
# ╠═bfb9d9d3-4560-4f49-86a0-01aafeaa37c8
# ╠═cf9df8e6-cf64-412c-a263-ba746f3ff706
# ╠═1eff07a7-c0a5-46f6-af57-9af2d45ddbd2
# ╠═5b9c1d07-d54b-4c3e-8962-f97793595d82
# ╠═3ce576c6-ab52-40ff-abad-b59521fd75a3
# ╠═52c4f0a4-c367-49b7-be89-5d038b5586aa
# ╠═fba6648b-b279-48cd-a956-0631eafc6899
# ╠═6d73d6bc-cb3d-488f-a47a-18aceba60dbf
# ╠═197fe84b-7314-4baf-ae79-4985b218db62
# ╠═96b32c90-2e8f-43d0-9407-3206f9ad9f0e
# ╠═a97405ed-6a50-4b97-9af2-5b8a074eabba
# ╠═47164d93-e263-405c-a071-e32270b83fd8
# ╠═29c34604-4b4b-4cf7-be40-ea8054ca173b
# ╠═1c390019-a2d4-4a2d-8db6-f0159a5ff82a
# ╠═c5a36cc9-5622-4abe-9dc4-5adbf0b999e9
# ╠═4daade11-9081-4c8c-a42a-7cda945c2611
# ╠═e29bf4c8-c686-4565-9c2a-33a1eb2b70c8
# ╠═fef63f90-b890-4a39-b577-aa420481b185
# ╠═6e9b4d7b-ff90-479b-aecf-a269f7716699
# ╠═3a0443be-f4f9-475c-b17a-2fe2ed32e7c4
