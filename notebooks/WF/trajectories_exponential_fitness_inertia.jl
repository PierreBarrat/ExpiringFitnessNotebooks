### A Pluto.jl notebook ###
# v0.19.22

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
datdir = "data_trajectories_exponential_fitness/"

# ╔═╡ e999242b-87b3-44a6-a264-9f1ffffc6101
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ ca412f4b-a14f-438a-8e09-7bc76c3467ff
βvals = map(x -> x["β"], values(files)) |> unique |> sort

# ╔═╡ d92d5654-2550-4c0f-b537-0b702cd2f860
αvals = map(x -> x["α"], values(files)) |> unique |> sort

# ╔═╡ 3cb788b1-ee3c-41c7-9663-72154a925aa3
fb = FrequencyBin(0.5, 0.05)

# ╔═╡ c60dc463-ab42-4943-b6e7-c3c0b0efc3bd


# ╔═╡ ead04896-0c4a-4e86-8521-668a522f1edc
1 /  0.005 / 0.6^2

# ╔═╡ 937244c8-7555-4f58-a715-1b7ce349a8d6
params = first(values(files))

# ╔═╡ b16e5600-d9f0-46f7-ab43-37421060b87b
# ╠═╡ skip_as_script = true
#=╠═╡
1+1
  ╠═╡ =#

# ╔═╡ f57bc80a-ab12-4911-bf98-3b3eefcb8ead
md"## Functions"

# ╔═╡ 6773d67d-2492-48b3-8c66-ee1199148606
function read_trajectories(file)
	return @chain begin
		file
		datdir*_
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories
	end
end

# ╔═╡ cffe4bb0-378d-4a0a-b51f-7fe9a3230958
always_below=false

# ╔═╡ ed3510c5-56f5-4404-9095-6d37cf6ad274
function _plot_trajectories(trajectories, fb)
	tmin = -600
	tmax = 2500
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
function plot_trajectories(params, fb)
	trajectories = read_trajectories(params["trajectory_file"])

	p = _plot_trajectories(trajectories, fb)
	sρ = round(params["s"] / params["ρ"]; sigdigits=2)
	plot!(
		p,
		title = "<β> = $(params["β"]), α = $(params["α"]), s/ρ = $(sρ), Δt = $(params["Δt"])"
	)
end

# ╔═╡ 46587230-b730-4ff7-be1a-ab186e70e951
plt_dict = map(values(files)) do x
	# @info x["trajectory_file"]
	file = x["trajectory_file"]
	(x["α"], x["β"], x["Δt"]) => plot_trajectories(x, fb)
end |> Dict

# ╔═╡ 235b0068-51b6-458a-8cbd-4ffe82c35083
p1 = let
	Δt = 10
	plot(
		[plt_dict[α, β, Δt] for α in αvals for β in βvals]...,
		layout = grid(length(αvals), length(βvals)),
		size = (length(βvals) * 600, length(αvals) * 600)
	)
end

# ╔═╡ 642da2d5-d196-4fcf-936c-a9585a7df465
savefig(p1, datdir * "inertia_panel_f$(fb.f).png")

# ╔═╡ Cell order:
# ╠═eadd1d72-d1f5-11ed-098a-11e1b44a86a9
# ╠═a55b7fe0-798b-4eca-b795-ee1342721d08
# ╠═e999242b-87b3-44a6-a264-9f1ffffc6101
# ╠═ca412f4b-a14f-438a-8e09-7bc76c3467ff
# ╠═d92d5654-2550-4c0f-b537-0b702cd2f860
# ╠═3cb788b1-ee3c-41c7-9663-72154a925aa3
# ╠═c60dc463-ab42-4943-b6e7-c3c0b0efc3bd
# ╠═46587230-b730-4ff7-be1a-ab186e70e951
# ╠═235b0068-51b6-458a-8cbd-4ffe82c35083
# ╠═ead04896-0c4a-4e86-8521-668a522f1edc
# ╠═642da2d5-d196-4fcf-936c-a9585a7df465
# ╠═937244c8-7555-4f58-a715-1b7ce349a8d6
# ╠═b16e5600-d9f0-46f7-ab43-37421060b87b
# ╟─f57bc80a-ab12-4911-bf98-3b3eefcb8ead
# ╠═6773d67d-2492-48b3-8c66-ee1199148606
# ╠═be559210-e29b-4369-b9f4-5232782b48df
# ╠═ed3510c5-56f5-4404-9095-6d37cf6ad274
# ╠═cffe4bb0-378d-4a0a-b51f-7fe9a3230958
