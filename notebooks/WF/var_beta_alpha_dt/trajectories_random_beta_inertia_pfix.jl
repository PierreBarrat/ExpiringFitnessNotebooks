### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ eadd1d72-d1f5-11ed-098a-11e1b44a86a9
begin
	using Revise
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using FrequencyTrajectories
	using JSON3
	using Measures
	using Plots
	using Random
	using StatsBase
end

# ╔═╡ a55b7fe0-798b-4eca-b795-ee1342721d08
datdir = "data_trajectories_random_beta.jl/"

# ╔═╡ e999242b-87b3-44a6-a264-9f1ffffc6101
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 3cb788b1-ee3c-41c7-9663-72154a925aa3
fb = FrequencyBin(0.4, 0.05)

# ╔═╡ f57bc80a-ab12-4911-bf98-3b3eefcb8ead
md"## Functions"

# ╔═╡ 6773d67d-2492-48b3-8c66-ee1199148606
function read_trajectories(file)
	return @chain begin
		file
		datdir * _
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories_from_dataframe
	end
end

# ╔═╡ c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
dat_all = let
	df = DataFrame(values(files))
	df.trajectories = map(f -> read_trajectories(f), df.trajectory_file)
	df
end;

# ╔═╡ 7bb49e6c-5ec9-4333-bc94-8c8d867ced84
dat = begin
	# dat_all[dat_all.α .== maximum(αvals), :]; # filter interesting part
	dat_all[dat_all.Δt .== 10, :];
end;

# ╔═╡ d92d5654-2550-4c0f-b537-0b702cd2f860
αvals = dat_all.α |> unique |> sort

# ╔═╡ b7a4f81a-f630-4b61-a8b2-70440fe9216f
Δtvals = dat_all.Δt |> unique |> sort

# ╔═╡ 318cdf36-8d40-4d27-961b-d2c4d8ba61b8
mβvals = dat_all.mβ |> unique |> sort

# ╔═╡ 5f687fbd-bbcd-4524-9eff-11f9493112c7
names(dat)

# ╔═╡ cffe4bb0-378d-4a0a-b51f-7fe9a3230958
always_below=true

# ╔═╡ ed3510c5-56f5-4404-9095-6d37cf6ad274
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

# ╔═╡ be559210-e29b-4369-b9f4-5232782b48df
function plot_trajectories(trajectories, params, fb)
	p = _plot_trajectories(trajectories, fb)
	# vline!([params.Ne], line=(:black, :dash), label="Ne=$(params.Ne)")
	plot!(
		p,
		title = "<β> = $(round(params.mβ, sigdigits=2)), α = $(params.α), Δt = $(params.Δt)"
	)
end

# ╔═╡ 46587230-b730-4ff7-be1a-ab186e70e951
plt_dict = map(eachrow(dat)) do r
	(r.α, r.mβ) => plot_trajectories(r.trajectories, r, fb)
end |> Dict

# ╔═╡ 235b0068-51b6-458a-8cbd-4ffe82c35083
p1 = let
	plot(
		[plt_dict[α, mβ] for α in αvals for mβ in mβvals]...,
		layout = grid(length(αvals), length(mβvals)),
		size = (length(mβvals) * 600, length(αvals) * 600),
		margin = 10mm,
		titlefontsize=18,
		legendfontsize=18,
		
	)
end

# ╔═╡ 642da2d5-d196-4fcf-936c-a9585a7df465
savefig(p1, datdir * "inertia_panel_f$(fb.f).png")

# ╔═╡ 54342047-a26e-4d2f-8f42-9f6bbb8f7816
function pfix_v_f(trajectories)
	fbs = [FrequencyBin(f, 0.05) for f in .1:.1:.9]
	pfix = 	map(fbs) do fb
		fixation_probability(T -> FT.inbin(T, fb; always_below), trajectories)
	end
	return (fbs, pfix)
end

# ╔═╡ 64bf3d22-666d-49af-9ba0-70a3858bc5fb
function plot_pfix(trajectories, params)
	fbs, pfix = pfix_v_f(trajectories)

	p = plot([fb.f for fb in fbs], pfix, label="", marker=:o)
	plot!([0,1], [0,1], line=(:black, :dash), label="")
	plot!(
		p,
		title = "<β> = $(round(params.mβ, sigdigits=2)), α = $(params.α), Δt = $(params.Δt)"
	)
end

# ╔═╡ 8e892a2e-c51d-405c-aa3a-3036bcf36a2c
plt_dict_pfix = map(eachrow(dat)) do r
	(r.α, r.mβ) => plot_pfix(r.trajectories, r)
end |> Dict

# ╔═╡ c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
p2 = let
	plot(
		[plt_dict_pfix[α, mβ] for α in αvals for mβ in mβvals]...,
		layout = grid(length(αvals), length(mβvals)),
		size = (length(mβvals) * 600, length(αvals) * 600)
	)
end

# ╔═╡ Cell order:
# ╠═eadd1d72-d1f5-11ed-098a-11e1b44a86a9
# ╠═a55b7fe0-798b-4eca-b795-ee1342721d08
# ╠═e999242b-87b3-44a6-a264-9f1ffffc6101
# ╠═c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
# ╠═7bb49e6c-5ec9-4333-bc94-8c8d867ced84
# ╠═d92d5654-2550-4c0f-b537-0b702cd2f860
# ╠═b7a4f81a-f630-4b61-a8b2-70440fe9216f
# ╠═318cdf36-8d40-4d27-961b-d2c4d8ba61b8
# ╠═3cb788b1-ee3c-41c7-9663-72154a925aa3
# ╠═46587230-b730-4ff7-be1a-ab186e70e951
# ╠═235b0068-51b6-458a-8cbd-4ffe82c35083
# ╠═642da2d5-d196-4fcf-936c-a9585a7df465
# ╠═8e892a2e-c51d-405c-aa3a-3036bcf36a2c
# ╠═c13f9cca-a534-4ecf-9b7a-0f8d4f825e55
# ╟─f57bc80a-ab12-4911-bf98-3b3eefcb8ead
# ╠═6773d67d-2492-48b3-8c66-ee1199148606
# ╠═5f687fbd-bbcd-4524-9eff-11f9493112c7
# ╠═be559210-e29b-4369-b9f4-5232782b48df
# ╠═ed3510c5-56f5-4404-9095-6d37cf6ad274
# ╠═64bf3d22-666d-49af-9ba0-70a3858bc5fb
# ╠═54342047-a26e-4d2f-8f42-9f6bbb8f7816
# ╠═cffe4bb0-378d-4a0a-b51f-7fe9a3230958
