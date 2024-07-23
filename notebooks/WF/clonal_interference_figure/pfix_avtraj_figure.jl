### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# ╔═╡ 078557cf-75b9-4c19-a6d0-7dbc5cdda594
begin
	using Pkg; Pkg.activate("../../../")
	using CSV
	using Chain
	using DataFrames
	using DataFramesMeta
	using Distributions
	using LaTeXStrings
	using Measures
	using QuadGK
	using StatsPlots
	using StatsBase
end

# ╔═╡ 8253380a-3dd2-11ef-050a-b736bacc969e


# ╔═╡ 88d61447-3cd0-448d-9cec-fe815bad464c
function pubfig(fnt_size=30; kwargs...)
    PLOTS_DEFAULTS = Dict(
        :markersize => 10,
        :linewidth => 5,
        :titlefontsize => fnt_size,
        :guidefontsize => fnt_size,
        :tickfontsize => fnt_size,
        :legendfontsize => fnt_size,
        :size => (1200,900),
        :gridlinewidth => 0.5,
        :gridalpha => 0.3,
        :framestyle => :box,
        :margin => 5mm,
        # :bottom_margin => 5mm,
        # :left_margin => 5mm,
    )

    for (k, x) in kwargs
        PLOTS_DEFAULTS[k] = x
    end

    return PLOTS_DEFAULTS
end

# ╔═╡ e24207be-c13e-4ef7-9fe5-7e56f9475646
begin
	plt_defaults = pubfig(24)
	Plots.default(; plt_defaults...)
end

# ╔═╡ 66cf7207-f348-46b4-b0b5-8d5118f3632a
cfs = "random"

# ╔═╡ 58f0eb7f-d658-4991-b832-3e4ce5b4307e
md"# Pfix"

# ╔═╡ 2fb56391-6cb2-4756-bd8d-aa0e3e6fd08d
md"# Av traj."

# ╔═╡ 2c86135b-c0c0-4814-9705-4792d9cd70f9
md"## Utils"

# ╔═╡ a73e555d-44d4-4b55-98c4-e6136b31a609
S_to_vec(S) = @chain split(S, r",| |\[|\]") filter!(!isempty, _) parse.(Float64, _)

# ╔═╡ 364c1696-c2e4-4de3-9f02-9c18bebeb163
data = let
	df = DataFrame(CSV.File("data_pfix_$(cfs).csv"))
	select!(
		df, Not(:pfix, :f), 
		[:f, :pfix] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:f, :pfix]
	)
	df.α = zeros(Float64, size(df, 1))
	df
end;

# ╔═╡ cdd7abd3-2dc0-4825-bcf8-4c363265135e
begin
	ρvals = data.ρ |> sort |> unique
	s0 = data.s[1]
end

# ╔═╡ 10b40c68-4c62-4f5e-9708-9cbd7a88c096
pal = palette(:bluesreds, length(ρvals))

# ╔═╡ 7c381b37-bd06-4f06-98b0-ef5753f66dd9
plot_pfix = let p = plot()
	for (k, ρ) in enumerate(ρvals)
		dat = @subset data :ρ .== ρ
		CI = round(ρ/s0, sigdigits=2)
		@df dat plot!(
			p, :f, :pfix;
			label = "ρ/s = $CI", 
			marker = :cross,
			linewidth = 5,
			markersize = 10,
			color = pal[k]
		)
	end
	plot!([0,1], [0,1], line=(:black, :dash), label="")
	p = plot!(
		xlabel = "frequency",
		frame=:box,
		ylabel = "P fixation",
		dpi = 300,
	)
	savefig("pfix_clonalinterference_regimes.png")
	p
end

# ╔═╡ 27dc6fa7-5eb2-493d-bc19-e74b584e010e
data_avtraj = let
	df = DataFrame(CSV.File("data_avtraj_$(cfs).csv"))
	select!(
		df, Not(:mt, :mf), 
		[:mt, :mf] => ByRow((x1, x2) -> (S_to_vec(x1), S_to_vec(x2))) => [:mt, :mf]
	)
	df.α = zeros(Float64, size(df, 1))
	df
end;

# ╔═╡ 3133fec3-5243-4317-955a-44f9a684a184
let	
	p = plot(
		xlabel = "time",
		ylabel = "frequency",
		frame=:box,
		title = "",
		# legend = (k == 1 ? :topleft : false),
		xlim = (-100, 1000),
		ylim = (-0.01, 1.01),
	)
		
	for (i, r) in enumerate(eachrow(data_avtraj))
		plot!(
			p, r.mt, r.mf;
			label="ρ/s=$(round(r.ρ/r.s; sigdigits=2))",
			linewidth = 2,
			color = pal[i],
		)
	end

	p
end

# ╔═╡ c07727bc-c72b-4c13-887a-930154f69d73
function estimate_pfix(x, s0, α)
	P(β) = (1-β)^(α/s0 - 1)
	eps = 1e-4
	Z = quadgk(β -> P(β), eps, 1-eps)[1]

	return x * quadgk(β -> P(β)/Z, eps, x)[1] + quadgk(β -> β*P(β)/Z, x, 1-eps)[1]
end

# ╔═╡ Cell order:
# ╠═8253380a-3dd2-11ef-050a-b736bacc969e
# ╠═078557cf-75b9-4c19-a6d0-7dbc5cdda594
# ╠═88d61447-3cd0-448d-9cec-fe815bad464c
# ╠═e24207be-c13e-4ef7-9fe5-7e56f9475646
# ╠═66cf7207-f348-46b4-b0b5-8d5118f3632a
# ╟─58f0eb7f-d658-4991-b832-3e4ce5b4307e
# ╠═364c1696-c2e4-4de3-9f02-9c18bebeb163
# ╠═cdd7abd3-2dc0-4825-bcf8-4c363265135e
# ╠═10b40c68-4c62-4f5e-9708-9cbd7a88c096
# ╠═7c381b37-bd06-4f06-98b0-ef5753f66dd9
# ╠═2fb56391-6cb2-4756-bd8d-aa0e3e6fd08d
# ╠═27dc6fa7-5eb2-493d-bc19-e74b584e010e
# ╠═3133fec3-5243-4317-955a-44f9a684a184
# ╠═2c86135b-c0c0-4814-9705-4792d9cd70f9
# ╠═a73e555d-44d4-4b55-98c4-e6136b31a609
# ╠═c07727bc-c72b-4c13-887a-930154f69d73
