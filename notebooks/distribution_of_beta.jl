### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 768fd5f6-ced3-11ed-18d3-fb9485be70ad
begin
	using Pkg; Pkg.activate("../")
	using Distributions
	using Plots
	using PlutoUI
	using StatsBase
end

# ╔═╡ 1c910626-d2f4-4896-9272-fffbbe0b0146
α = .01

# ╔═╡ 3280a381-96e3-4451-96d5-09949e405a6a
s0 = 0.005

# ╔═╡ 228e85cc-269e-44db-aa60-84210d2b4537
β(s, α) = 1 - exp(-s/α)

# ╔═╡ 55129e4b-ae94-4cb9-8323-bec17a8b933d
fitness_distribution = let
	Exponential(s0)
	# Dirac(1)
end

# ╔═╡ f98f4eb5-b570-4bac-81c7-82400824537e
begin
	s_values = rand(fitness_distribution, 10_000)
	β_values = map(s -> β(s, α), s_values)
end

# ╔═╡ 6e5ceae4-8229-45f6-8d9c-a2568f146231
let
	bins = 0:0.01:1
	bin_centers = (bins[2:end] + bins[1:end-1])/2
	h = fit(Histogram, β_values, bins)
	plot(bin_centers, h.weights, xlabel="β")
end

# ╔═╡ 6b44f2ba-7a0e-4cf8-b0e3-d9bed51dcdda
md"### Moments"

# ╔═╡ 9308dc8e-7b78-46ba-8fbf-ca91a36a1ed9
mean(β_values.^2)

# ╔═╡ 6b2a5fb6-fea9-4da1-b80e-1c64456de16c
s0 / (s0 + α) * (2*s0)/ (2*s0 + α)

# ╔═╡ c8114c7c-44f1-42fd-ac97-29a98422be25
mean(β_values)

# ╔═╡ 56b773df-a208-436d-9715-04f43eee3984
s0 / (s0 + α)

# ╔═╡ 689d852a-09b1-4950-87d1-34264a576ed6
md"## parametrizing with first and second moment"

# ╔═╡ 57cdf9ab-9fbd-4e64-9cc3-6f056479d829
md"""
The Beta distribution is traditionally parametrized by two numbers $a, b>0$, with first and second moments being

$$\langle\beta\rangle = m = \frac{a}{a+b}\;\;\text{and}\;\;\langle\beta^2\rangle\frac{a}{a+b}\frac{a+1}{a+b+1}.$$

For our purpose, it is easier to parametrize it using the first and second moments. This can be done using the following identities: 

$$\begin{align}
b &= \frac{\langle\beta\rangle}{\langle\beta^2\rangle - \langle\beta\rangle^2}(1-\langle\beta\rangle),\\
a &= \frac{\langle\beta\rangle}{1-\langle\beta\rangle}b.
\end{align}$$

Note that there are the following bounds for the moments as parameters: 

$$0 \leq \langle\beta^2\rangle \leq \langle\beta\rangle \leq \sqrt{\langle\beta^2\rangle} \leq 1.$$

Or stated in another way: 

$$\begin{align}
\langle\beta\rangle^2 &\leq \langle\beta^2 \rangle \leq \langle\beta\rangle \\
\langle\beta^2\rangle &\leq \langle \beta \rangle \leq \sqrt{\langle\beta^2\rangle}
\end{align}$$
"""

# ╔═╡ e71073d8-e76d-45d1-aa67-9ffdad6c98e3
function get_β_distribution(mβ, β2)
	b = (mβ - β2)/(β2 - mβ^2)*(1-mβ)
	a = mβ/(1-mβ)*b
	return Beta(a,b)
end

# ╔═╡ 6bfb1a84-3b30-47aa-9fa7-ea5218142496
begin
	Ne = 100
	ρ = 0.1
	β2 = 1/ρ/Ne
end

# ╔═╡ 65fbcfe1-7b31-4c14-a129-d5290cae357a
mβvals = collect(range(β2, sqrt(β2), length=7))[2:end-1]

# ╔═╡ 4b5372f1-67ab-417f-a57f-9df130e4f4b0
plts = map(mβvals) do mβ
	B = get_β_distribution(mβ, β2)
	x = 0.:0.01:1
	plot(
		x, pdf.(B, x), label="",
		title = "<β> = $(round(mβ, sigdigits=2)), <β^2>=$(β2)"
	)
end

# ╔═╡ 38d32a2b-1c35-4a20-9bb2-ad302fcfa192
let
	p = plot(plts..., layout=(1,5), size = (3000, 600))
	savefig("Beta_distribution_fixed_beta2.png")
	p
end

# ╔═╡ 5234c474-80eb-45fb-a3f6-a14b2f7753bd
@bind mβ Slider(0.05:0.05:0.95, show_value=true, default = 0.5)

# ╔═╡ 340b7282-bbbd-4649-8eea-c58b7d850f16
let
	β2 = mβ^2 * 1.01
	B = get_β_distribution(mβ, β2)
	x = 0.:0.001:1
	plot(x, pdf.(B, x), label="β2=$(β2)",)

	β2 = round(mβ^2 + (mβ - mβ^2)/3, sigdigits=3)
	B = get_β_distribution(mβ, β2)
	x = 0.:0.001:1
	plot!(x, pdf.(B, x), label="β2=$(β2)",)
	
	# β2 = mβ * 0.9
	# B = get_β_distribution(mβ, β2)
	# x = 0.:0.01:1
	# plot!(x, pdf.(B, x), label="β2=$(β2)",)

	plot!(title = "<β> = $(round(mβ, sigdigits=2))")
end

# ╔═╡ 401095ef-af81-40f1-b914-8c9a048d45a7
md"## Mass above $1-e^{-3\rho/\alpha}$"

# ╔═╡ 6ab4ebdc-3204-46c0-b357-cdb148619d4c
0.2:0.15:0.65 |> collect 

# ╔═╡ 7e7c9997-7f2a-415b-8ed5-fb3719d15190
500/0.05

# ╔═╡ d6e7301a-cf3f-4b49-ac25-077e5b9eaf5a
plts_2 = let
	ρ = .05
	α = .5
	thr = 1 - exp(-3*ρ/α)
	
	mβvals = 0.2:0.15:0.65 |> collect
	β2vals = map(mβvals) do m
		r = range(m^2*1.01, m*0.99, length=8)
		L = div(length(r), 2)
		r[1:L]
	end
	x = length(mβvals)
	y = length(first(β2vals))
	

	plts = Matrix{Any}(undef, y, x)
	for (i, m) in enumerate(mβvals), (j,v) in enumerate(β2vals[i])
		B = get_β_distribution(m, v)
		below = round(cdf(B, thr), sigdigits=2)
		local x = 0.:0.001:1
		p = plot(
			x, pdf.(B, x), label="$below",
			title = "<β> = $(round(m, sigdigits=2)), <β^2>=$(round(v, sigdigits=2))",
		)
		vline!([thr], line=(:black, :dash), label="")
		plts[j,i] = p
	end
	
	plot(plts..., layout=(x, y), size = (1200, 1200))
	
end

# ╔═╡ e858d9ff-63e4-4913-9393-3c050ee74889
let
	ρ = .05
	α = .5
	thr = 1 - exp(-3*ρ/α)
end

# ╔═╡ 02046fbd-eca7-4ea5-9426-0b7129830608
md"# Tests"

# ╔═╡ bc9620d5-a13d-4158-9dd5-79b889e2c541
let
	mβ = 0.7
	β2 = 0.525
	B = get_β_distribution(mβ, β2)
	x = rand(B, 10_000)

	plot(sort(x), 1:length(x), title = mean(x))
	plot(0:0.01:1, pdf.(B, 0:0.01:1))
end

# ╔═╡ Cell order:
# ╠═768fd5f6-ced3-11ed-18d3-fb9485be70ad
# ╠═1c910626-d2f4-4896-9272-fffbbe0b0146
# ╠═3280a381-96e3-4451-96d5-09949e405a6a
# ╠═228e85cc-269e-44db-aa60-84210d2b4537
# ╠═55129e4b-ae94-4cb9-8323-bec17a8b933d
# ╠═f98f4eb5-b570-4bac-81c7-82400824537e
# ╠═6e5ceae4-8229-45f6-8d9c-a2568f146231
# ╟─6b44f2ba-7a0e-4cf8-b0e3-d9bed51dcdda
# ╠═9308dc8e-7b78-46ba-8fbf-ca91a36a1ed9
# ╠═6b2a5fb6-fea9-4da1-b80e-1c64456de16c
# ╠═c8114c7c-44f1-42fd-ac97-29a98422be25
# ╠═56b773df-a208-436d-9715-04f43eee3984
# ╠═689d852a-09b1-4950-87d1-34264a576ed6
# ╟─57cdf9ab-9fbd-4e64-9cc3-6f056479d829
# ╠═e71073d8-e76d-45d1-aa67-9ffdad6c98e3
# ╠═6bfb1a84-3b30-47aa-9fa7-ea5218142496
# ╠═65fbcfe1-7b31-4c14-a129-d5290cae357a
# ╠═4b5372f1-67ab-417f-a57f-9df130e4f4b0
# ╠═38d32a2b-1c35-4a20-9bb2-ad302fcfa192
# ╠═5234c474-80eb-45fb-a3f6-a14b2f7753bd
# ╠═340b7282-bbbd-4649-8eea-c58b7d850f16
# ╟─401095ef-af81-40f1-b914-8c9a048d45a7
# ╠═6ab4ebdc-3204-46c0-b357-cdb148619d4c
# ╠═7e7c9997-7f2a-415b-8ed5-fb3719d15190
# ╠═d6e7301a-cf3f-4b49-ac25-077e5b9eaf5a
# ╠═e858d9ff-63e4-4913-9393-3c050ee74889
# ╠═02046fbd-eca7-4ea5-9426-0b7129830608
# ╠═bc9620d5-a13d-4158-9dd5-79b889e2c541
