### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 7d7eba62-d2f2-11ed-3965-21bef91270c2
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

# ╔═╡ db84f3b1-28d4-4bbe-941d-f799117536b9
datdir = "data_length_trajectories_exponential_fitness/"

# ╔═╡ 92d2b993-ccbc-4ac6-a735-f16eb836571a
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 0f2f33d4-e8e0-42fc-971c-4b9aa36c2699
βvals = map(x -> x["β"], values(files)) |> unique |> sort

# ╔═╡ 926d0d74-6ce0-402f-a236-372d4a1ef355
ρvals = map(x -> x["ρ"], values(files)) |> unique |> sort

# ╔═╡ 68659fb3-1637-4cf4-b2ee-6416b9465586
begin
	α = first(values(files))["α"]
	L = first(values(files))["L"]
	N = first(values(files))["N"]
	Δt = first(values(files))["Δt"]
end

# ╔═╡ 65364020-7680-4019-bb4c-d0651f4a041b
always_below = true

# ╔═╡ 4bef6c82-8207-4093-9910-a1ad9a6d8d4f
fb = FrequencyBin(0.3, 0.05)

# ╔═╡ fb2ba0f9-40fc-43df-86d2-a1ab1f6f6417
md"## Functions"

# ╔═╡ 191a9bba-6c43-4d97-a875-f3441ba86e41
function βeff(s, ρ, α, x0) 
	β = 1 - exp(-s/α)
	x0*β*exp(s/ρ) / (x0 * (exp(s/ρ) - 1) + β)
end

# ╔═╡ 4cab3539-3978-4117-bc38-0f039217944d
function read_trajectories(file)
	return @chain begin
		file
		datdir*_
		CSV.read(_, DataFrame) 
		FrequencyTrajectories.trajectories
	end
end

# ╔═╡ af441dd6-4274-4d09-ac75-bc6b0bdfe7dc
function _hist_trajectory_length(trajectories; bins = :auto)
	
	local p = plot(
		xlabel = "Length (gen.)",
		size = (900, 900),
		frame = :box
	)
	X = map(FrequencyTrajectories.duration, trajectories)
	histogram!(X; bins)
	
	return p
end

# ╔═╡ 6993f547-011d-4a83-bd5e-3c7273ba1d42
function hist_trajectory_length(params; fb=nothing)
	file = params["trajectory_file"]
	Ne = params["Ne"]
	
	trajectories = @chain begin
		read_trajectories(file)
		isnothing(fb) ? _ : filter(T -> FT.inbin(T, fb), _)
		filter(T -> T.final_state != :missing, _)
	end

	bins = range(0, 10*Ne, 25)
	
	p = _hist_trajectory_length(trajectories; bins)
	vline!([Ne], line=(5, 0.6, :black, :dash))
	plot!(
		p,
		title = "ρ = $(params["ρ"]), <β> = $(params["β"]), Ne = $(round(Ne, sigdigits=3)), s/ρ = $(round(params["s"]/params["ρ"]; sigdigits=2))"
	)

	Ne_eff = 1/params["ρ"]/βeff(params["s"], params["ρ"], params["α"], 0.02)^2
	# vline!([Ne_eff], line=(5, 0.6, :black))

	return p
end

# ╔═╡ 6b92286f-d1ac-44e5-a815-0a6b4219a0b6
plt_dict = map(values(files)) do p
	(p["ρ"], p["β"]) => hist_trajectory_length(p; fb)
end |> Dict

# ╔═╡ 9c0b0021-fb3a-4e7d-ba8d-675edc59afbd
p = plot(
	[plt_dict[ρ, β] for ρ in ρvals for β in βvals]...,
	layout = grid(length(ρvals), length(βvals)),
	size = (length(βvals) * 600, length(ρvals) * 600)
)

# ╔═╡ 5247077d-aae3-4295-9845-0791bda6d0df
savefig(p, datdir*"length_histograms.png")

# ╔═╡ 6b071874-6f2a-45da-ad91-a5e696b9ec23
md"## Tests"

# ╔═╡ d56a021c-cca6-412d-a4b6-3495bb70967f
let
	s = .05
	α = .1
	ρ = .1

	β = 1-exp(-s/α)

	x0 = 1e-5
	tvals = range(0, 12/s, 1000)
	x = map(tvals) do t 
		β * x0 * exp(s*t) / (x0*(exp(s*t)-1) + β)
	end
	plot(tvals, x, label="")
	hline!([β], line=(:black, :dash), label="β")
	vline!([1/ρ], label = "1/ρ")
	hline!([β * x0 * exp(s/ρ) / (x0*(exp(s/ρ)-1) + β)], label="effective β")

	plot!(legend = :outerright)
end
	

# ╔═╡ 03255ba2-868f-43aa-b637-46d5f718636c
let
	β = .3
	α = .1
	s = -α * log(1-β)

	ρvals = s * 10 .^ range(-2, 1, length=100)
	
	p = plot(
		xlabel = "ρ/s",
		title = "effective β",
		xscale=:log10,
	)

	x0 = .02
	plot!(ρvals/s, βeff.(s, ρvals, α, x0), label="x0 = 0.02", )

	x0 = 1e-5
	plot!(ρvals/s, βeff.(s, ρvals, α, x0), label="x0 = 1e-5", )
end

# ╔═╡ fdd7cea3-6a31-4dfb-a0ca-b9924c59db23
1e-5 * 0.3 * exp(10) / (1e-5*exp(10) + 0.3)

# ╔═╡ 5317d975-a301-42e7-9a63-7aa32b4ac4fe
let
	β = .6
	α = .1

	s = -α * log(1-β)

	ρvals = s * 10 .^ range(-2, 1, length=100)
	p = plot(
		xlabel = "ρ/s",
		title = "effective β",
		xscale=:log10,
	)

	x0 = .02
	plot!(ρvals/s, βeff.(s, ρvals, α, x0), label="x0 = 0.02", )

	x0 = 1e-5
	plot!(ρvals/s, βeff.(s, ρvals, α, x0), label="x0 = 1e-5", )
end

# ╔═╡ Cell order:
# ╠═7d7eba62-d2f2-11ed-3965-21bef91270c2
# ╠═db84f3b1-28d4-4bbe-941d-f799117536b9
# ╠═92d2b993-ccbc-4ac6-a735-f16eb836571a
# ╠═0f2f33d4-e8e0-42fc-971c-4b9aa36c2699
# ╠═926d0d74-6ce0-402f-a236-372d4a1ef355
# ╠═68659fb3-1637-4cf4-b2ee-6416b9465586
# ╠═65364020-7680-4019-bb4c-d0651f4a041b
# ╠═4bef6c82-8207-4093-9910-a1ad9a6d8d4f
# ╠═6b92286f-d1ac-44e5-a815-0a6b4219a0b6
# ╠═9c0b0021-fb3a-4e7d-ba8d-675edc59afbd
# ╠═5247077d-aae3-4295-9845-0791bda6d0df
# ╟─fb2ba0f9-40fc-43df-86d2-a1ab1f6f6417
# ╠═191a9bba-6c43-4d97-a875-f3441ba86e41
# ╠═4cab3539-3978-4117-bc38-0f039217944d
# ╠═6993f547-011d-4a83-bd5e-3c7273ba1d42
# ╠═af441dd6-4274-4d09-ac75-bc6b0bdfe7dc
# ╠═6b071874-6f2a-45da-ad91-a5e696b9ec23
# ╠═d56a021c-cca6-412d-a4b6-3495bb70967f
# ╠═03255ba2-868f-43aa-b637-46d5f718636c
# ╠═fdd7cea3-6a31-4dfb-a0ca-b9924c59db23
# ╠═5317d975-a301-42e7-9a63-7aa32b4ac4fe
