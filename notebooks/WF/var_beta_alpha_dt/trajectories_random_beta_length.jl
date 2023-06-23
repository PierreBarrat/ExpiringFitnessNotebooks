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

# ╔═╡ 35e00ed3-5acf-41f6-908f-1190768f6843
function read_trajectories(file)
	CSV.read(file, DataFrame) |> FrequencyTrajectories.trajectories_from_dataframe
end

# ╔═╡ 10987b07-70e0-4a19-9173-83cf487c08ff
dat_all = let
	df = DataFrame(values(files))
	df.trajectories = map(f -> read_trajectories(datdir * f), df.trajectory_file)
	df
end;

# ╔═╡ 2b33fcb7-9270-42a6-8c20-74178d8acfcb
dat = dat_all[dat_all.Δt .== 1, :];

# ╔═╡ 5d97c89c-1eb5-4f57-b3aa-012328964944
names(dat)

# ╔═╡ 3fc64542-ced7-4c92-9c60-6e88421e7374
αvals = dat_all.α |> unique |> sort

# ╔═╡ 8146331b-01a7-49b0-abce-b395824b4a7a
mβvals = dat_all.mβ |> unique |> sort

# ╔═╡ 7ba85594-bf3c-4179-9cb3-78ee4ea1d9a9
function length_histogram(trajectories, params)
	X = map(FT.duration, trajectories)
	bins = range(0, maximum(X)*1.1, length=20)
	h = fit(Histogram, map(FT.duration, trajectories), bins)

	bin_centers = (bins[2:end] + bins[1:end-1])/2
	plot(bin_centers, h.weights / sum(h.weights), label="", linewidth=3)
	
	Ne = round(1/params.ρ/params.β2, sigdigits = 2)
	vline!([Ne], line=(:black, :dash), linewidth=3, label="Ne=$Ne")
	vline!([mean(X)], label="", color=1, line=(:dash), linewidth=3)

	plot!(title = "<β>=$(round(params.mβ, sigdigits=2)), α=$(params.α)", xlabel = "time")

	plot!(ylim = (-0.05, 0.5))
end

# ╔═╡ f566f619-e44d-47d8-9506-9df7aa067022
plt_dict = map(eachrow(dat)) do r
	(r.α, r.mβ) => length_histogram(r.trajectories, r)
end |> Dict

# ╔═╡ 6c9414ed-bfd7-4afb-a6a9-7ca373e81fe3
p1 = let
	plot(
		[plt_dict[α, mβ] for α in αvals for mβ in mβvals]...,
		layout = grid(length(αvals), length(mβvals)),
		size = (length(mβvals) * 600, length(αvals) * 600),
		margin=10mm,
		titlefontsize=20,
		legendfontsize=20,
		tickfontsize=16
	)
end

# ╔═╡ Cell order:
# ╠═5127cf4e-eb40-11ed-0a15-53a1acf43aef
# ╠═5ffd9c4e-4807-43aa-943a-4c8c3f9b022a
# ╠═2713455c-df95-49c7-925a-9d92b91719a5
# ╠═10987b07-70e0-4a19-9173-83cf487c08ff
# ╠═2b33fcb7-9270-42a6-8c20-74178d8acfcb
# ╠═5d97c89c-1eb5-4f57-b3aa-012328964944
# ╠═3fc64542-ced7-4c92-9c60-6e88421e7374
# ╠═8146331b-01a7-49b0-abce-b395824b4a7a
# ╠═f566f619-e44d-47d8-9506-9df7aa067022
# ╠═6c9414ed-bfd7-4afb-a6a9-7ca373e81fe3
# ╠═9ac16620-e45c-4218-9d7f-c7cca89e66e1
# ╠═35e00ed3-5acf-41f6-908f-1190768f6843
# ╠═7ba85594-bf3c-4179-9cb3-78ee4ea1d9a9
