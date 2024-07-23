### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 2168f2c6-4905-11ef-35d0-515b0eba30d2
begin
	using Pkg; Pkg.activate("../../")
	using Distributions
	using Plots
	using Measures
	using StatsBase
end

# ╔═╡ 5ac8b560-da86-4fc8-9b11-096e72d96db0
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ b817799b-248c-43c6-b4b2-09207fc136d1
let
	font_size = 22
	Plots.default(; pubfig(font_size)...)
end

# ╔═╡ 90c3867c-8e76-40f9-83f9-5f6b894f518a
md"# Exponential distribution of $(1-f)$"

# ╔═╡ 2fd52b83-df6f-4ce5-905c-5caec3f27037
function Pβ_exp(β, λ, μ)
	return μ/λ / (μ/λ*β + (1-β))^2
end

# ╔═╡ da956c27-782b-4974-9e32-41fa025b3293
function mβ_exp(λ, μ)
	xvals = range(0, 1, length = 10_000)
	mβ = sum(x -> Pβ_exp(x, λ, μ)*x, xvals)/length(xvals)
	vβ = sum(x -> Pβ_exp(x, λ, μ)*x^2, xvals)/length(xvals) - mβ^2
	return mβ, vβ
end

# ╔═╡ 6e5cdc2c-d146-4b3c-adb5-6962eabbb50f
begin
	λ = .05
	μvals = [.01, .025, .05, .1]
end

# ╔═╡ fd30fc51-c3cb-4446-bbbe-a095cf5e5348
mβ_exp(λ, .001)

# ╔═╡ d6ded531-07f1-4bfc-9104-7d670578ff57
plt_pβ_exp = let p = plot()
	pal = palette(:redsblues, length(μvals))
	
	βvals = range(0, 1, length = 1000)
	for (i, μ) in enumerate(μvals) 
		plot!(
			βvals, map(x -> Pβ_exp(x, λ, μ), βvals);
			label = "μ/λ=$(round(μ/λ; sigdigits=2))", color = pal[i], 
		)
	end
	plot!(
		xlabel = "β",
		ylabel = "P(β)",
	)
end

# ╔═╡ a7004b28-a097-46d6-9141-430c9fe45396
logrange(a,b;length=2) = exp.(range(log(a), log(b); length))

# ╔═╡ 5d5bdb36-c0c5-4ea4-940a-eb68e6b63530
plt_mβ_exp = let p = plot()
	μvals = logrange(λ/100, λ*100, length = 1000)
	X = map(μ -> mβ_exp(λ, μ), μvals)
	mβ, vβ = ([x[1] for x in X], [x[2] for x in X])
	plot!(μvals./λ, mβ; label = "<β>")
	plot!(μvals./λ, sqrt.(vβ); label = "std(β)")
	plot!(
		xlabel = "μ/λ",
		xscale = :log10,
		ylim = (0,1),
		xticks = [1e-2, 1e-1, 1e0, 1e1]
,	)
end

# ╔═╡ b0c96f79-85b9-4fae-aac1-0e7a2bb5726b
let
	p = plot(
		plt_pβ_exp, plt_mβ_exp; 
		layout = grid(1,2), size = (1200, 600), dpi = 300, 
		bottom_margin = 10mm, left_margin = 10mm,
	)
	savefig("../../Figures/beta_distribution_v_K.png")
	p
end

# ╔═╡ Cell order:
# ╠═2168f2c6-4905-11ef-35d0-515b0eba30d2
# ╠═5ac8b560-da86-4fc8-9b11-096e72d96db0
# ╠═b817799b-248c-43c6-b4b2-09207fc136d1
# ╟─90c3867c-8e76-40f9-83f9-5f6b894f518a
# ╠═2fd52b83-df6f-4ce5-905c-5caec3f27037
# ╠═da956c27-782b-4974-9e32-41fa025b3293
# ╠═6e5cdc2c-d146-4b3c-adb5-6962eabbb50f
# ╠═fd30fc51-c3cb-4446-bbbe-a095cf5e5348
# ╟─d6ded531-07f1-4bfc-9104-7d670578ff57
# ╟─5d5bdb36-c0c5-4ea4-940a-eb68e6b63530
# ╠═b0c96f79-85b9-4fae-aac1-0e7a2bb5726b
# ╠═a7004b28-a097-46d6-9141-430c9fe45396
