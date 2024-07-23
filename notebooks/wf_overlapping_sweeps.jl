### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 205722a2-159b-11ee-1d99-0106c03d3ad0
begin
	using Pkg; Pkg.activate("..")
	using Chain
	using Plots
	using StatsBase
	using WrightFisher
end

# ╔═╡ 73a7eeba-8a53-46e5-9dd7-181f93412572
begin
	N = 100_000_000
	L = 200
	Δt = 10
	
	α = 0.1
	β = 0.4
	s = -log(1-β)*α
end

# ╔═╡ 4709dca9-9739-4957-8bdf-e18700b225e6
function _simulate(ρ)
	evtime = max(15/s, 3/ρ)
	L = 2*round(Int, evtime * ρ + 10) |> Int
	ϕ = WF.ExpiringFitness(L, zeros(Float64, L), α)
	pop = WF.Pop(ϕ; N, L, μ=0)

	cb = (
		f1 = WF.frequencies,
	)
	
	cb_vals, switch_times = WF.Tools.evolve_sample!(
		pop, evtime, Δt, cb;
		fitness_distribution = s/2, 
		switchgen = 1/ρ, 
		change_init_field = false, 
		change_field_time = :periodic,
	)
	# @info "# switches:" length(switch_times)
	# @info "L: " L
	return cb_vals, switch_times
end

# ╔═╡ 594f3fff-b398-47c3-9edc-b99f657704a5
function simulate(ρ; Nrep = 100)
	# average trajectory of the first partial sweep 
	f1 = mean(1:Nrep) do _
		cb_vals, switch_times = _simulate(ρ)
		i_first = switch_times[1].pos
		map(x -> x.f1[2*i_first], cb_vals)
	end
	# just to get time values
	cb_vals, _ = _simulate(ρ)
	tvals = map(x -> x.t, cb_vals)
	
	return tvals, f1
end

# ╔═╡ 802122e8-8e13-4ac2-a8f5-94e1e4298ede


# ╔═╡ 76d49496-e85d-41ce-84e5-e154623b8aa4
k = 3

# ╔═╡ 904df873-9589-4881-9599-6b72294c8b90
β * (1-β)^k

# ╔═╡ c2c47bef-38a9-47e2-85de-1d00324d3fe2
ρ = k*s

# ╔═╡ 9f128504-e947-4050-844b-a19a96183431
βe = (1 - exp(-(1+k)*s/α))/(1+k)

# ╔═╡ e6ce0eec-8253-417d-8492-66b33ee525ec
1/ρ

# ╔═╡ ac271006-1636-47af-af85-5e07c3db1e32
logrange(x, y, L) = exp.(range(log(x), log(y), length=L))

# ╔═╡ af33976f-fe89-425b-9974-0a10712f0d27
ρ_values = logrange(0.01, 30, 15) * s

# ╔═╡ 02ed5447-f6d7-47c3-855e-0ede01f6b9be
map(ρ_values) do ρ
	evtime = max(15/s, 3/ρ)
	round(Int, evtime * ρ + 10) |> Int
end

# ╔═╡ ae55a216-2052-43b3-80b6-1093b11ce08d
final_values = map(ρ_values) do ρ
	_, f = simulate(ρ, Nrep = 2500)
	mean(f[end-10:end])
end

# ╔═╡ 23b2eb77-d30d-49bc-a610-55f6f374f2d2
let
	p = plot(
		xlabel = "ρ (scaled to s)",
		ylabel = "final frequency"
	)
	plot!(1 ./(ρ_values/s), final_values, label="", xscale=:log10)

	hline!([β], line = (:black, :dash))

	# kvals = map(ρ -> ρ/s, ρ_values)
	# plot!(ρ_values/s, β*(1-β).^kvals)
	# plot!(ρ_values/s, map(k -> (1 - exp(-(1+k)*s/α))/(1+k), kvals))
	# plot!(
	# 	ρ_values/s, 
	# 	map(ρ_values) do ρ
	# 		k = ρ/s
	# 		βe = (1 - exp(-(1+k)*s/α))/(1+k)
	# 		T = log(β/0.02)/s
	# 		(1-exp(-ρ*T)) * βe + exp(-ρ*T) * β
	# 	end
	# )
end

# ╔═╡ 6a8d5437-f79e-486f-bb17-ff06a0544680
tvals, f = simulate(ρ_values[1]; Nrep = 10)

# ╔═╡ 8dff80cc-6ec1-4e57-9211-be0f1edcf419
let
	plot(tvals*s, f, label="")
	plot!(ylim = (-0.025,1.025))
	# hline!([β], label="β")
	hline!([β*(1-β)^k], label="β(1-β)^k")
	hline!([βe], label="βe")
end

# ╔═╡ Cell order:
# ╠═205722a2-159b-11ee-1d99-0106c03d3ad0
# ╠═73a7eeba-8a53-46e5-9dd7-181f93412572
# ╠═904df873-9589-4881-9599-6b72294c8b90
# ╠═4709dca9-9739-4957-8bdf-e18700b225e6
# ╠═594f3fff-b398-47c3-9edc-b99f657704a5
# ╠═af33976f-fe89-425b-9974-0a10712f0d27
# ╠═02ed5447-f6d7-47c3-855e-0ede01f6b9be
# ╠═ae55a216-2052-43b3-80b6-1093b11ce08d
# ╠═23b2eb77-d30d-49bc-a610-55f6f374f2d2
# ╠═c2c47bef-38a9-47e2-85de-1d00324d3fe2
# ╠═802122e8-8e13-4ac2-a8f5-94e1e4298ede
# ╠═6a8d5437-f79e-486f-bb17-ff06a0544680
# ╠═76d49496-e85d-41ce-84e5-e154623b8aa4
# ╠═9f128504-e947-4050-844b-a19a96183431
# ╠═8dff80cc-6ec1-4e57-9211-be0f1edcf419
# ╠═e6ce0eec-8253-417d-8492-66b33ee525ec
# ╠═ac271006-1636-47af-af85-5e07c3db1e32
