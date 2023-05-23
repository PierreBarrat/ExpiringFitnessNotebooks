### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ c7e4d226-f881-11ed-0717-0f30a6dc429a
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using CSV
	using DataFrames
	using DelimitedFiles
	using KernelDensitySJ
	using JSON3
	using Measures
	using Plots
	using StatsBase
end

# ╔═╡ 385ae3ab-f8d3-4d8d-9cca-43730702083c
datdir = "data_trajectories_var_rho.jl/"

# ╔═╡ 6a7e0c08-f9e0-475b-80f7-ebc2a69294a6
files = open(datdir * "files.json", "r") do io
	JSON3.read(io, Dict)
end

# ╔═╡ 167cdfe9-fa90-4982-92f0-430906199130
ρvals = map(x -> x["ρ"], values(files)) |> unique |> sort

# ╔═╡ b2ccf546-18aa-49a3-8499-ac21d4cc54ce
mβ_vals = map(x -> x["mβ"], values(files)) |> unique |> sort

# ╔═╡ 37d544e6-e6a5-4e08-852c-d67820774fe5
α = map(x -> x["α"], values(files)) |> unique |> first

# ╔═╡ 13826c8a-1770-4f6c-a5a3-6b0163fafa82
md"""
## Time for a partial sweep

The partial sweep with fitness effect $s$ starting at frequency $x_0$ and ending at $\beta$ can be approximated by the function

$$x(t) = \frac{\beta x_0 e^{st}}{x_0 e^{st} + \beta - x_0}.$$

We have $x(0) = x_0$ and $x(\infty) = \beta$. 
The time needed for the sweep to reach $\nu\beta$ with $\nu < 1$ is  

$$\begin{align}
\frac{\beta x_0 e^{st}}{x_0 e^{st} + \beta - x_0} &= \nu\beta,\\
x_0e^{st} (1-\nu) &= \nu(\beta-x_0),\\
t &= s^{-1}\log\frac{\nu}{1-\nu}\frac{\beta-x_0}{x_0}.
\end{align}$$

Here I use $x_0 = 0.02$ and we can use $\nu=0.9$, meaning

$$t \simeq s^{-1} (6 + \log(\beta - x_0)).$$
"""

# ╔═╡ fa04fcc4-081f-467f-b2ac-b890be897c83
begin
	ν = 0.9
	x0 = 0.02
end

# ╔═╡ b8dea52e-be5c-47a3-9bc4-203fa360a035
md"""
Reading info about sweeps, and computing $\beta$, the time to the next sweep and the expected time for the actual sweep to complete. 
"""

# ╔═╡ 60700a60-d523-4429-90f5-62ea201c3dfd
md"""
Below, the time needed for a fraction $\nu$ of the partial sweep to be finished (using the fitness of the mutation) vs. the observed time before the next sweep. 
"""

# ╔═╡ 2f25a0eb-96a7-4e87-acd0-cdc4c10e1b05
md"## Change in frequency of the sweeping allele"

# ╔═╡ 8a929898-b50e-43d3-bf14-bf8e72ca8ee2
md"## Change in frequency of the other alleles"

# ╔═╡ b58c5cab-b329-4207-9f05-4e4ce0e9e98f
md"""
Assume there is a partial sweep of amplitude $\beta$ at position $i$ at time $t$. The question we ask is: what is the change in frequency $\Delta f$ of the alleles at positions $j\neq i$ after some time $\Delta t$. Consider $j \neq i$ with an allele frequency $x_0$. We only consider the caes $x_0 \in [0.05, 0.95]$. If $\Delta t$ is long enough, we expect  two cases: 
-  $\Delta f_j^{observed} > 0$, in which case we should expect $\Delta f^{theory} = \beta (1-x_0)$; the scaled change in frequency is then 
$$\Delta f_j^{scaled} = \frac{\Delta f_j^{observed}}{\beta(1-x_0)};$$
-  $\Delta f_j^{observed} < 0$, in which case we have with the same logic
$$\Delta f_j^{scaled} = -\frac{\Delta f_j^{observed}}{\beta x_0}.$$

Note that we can choose two values for $\Delta t$: the time when the next sweep starts, or the time it would take the partial sweep to almost complete if unperturbed. 

For each partial sweep with size $\beta$, we plot $\Delta f_j ^{scaled}$ for all $j$, and expect the values to land on $1$ if the sweeps do not overlap. 
"""

# ╔═╡ 1fe37cb1-b70a-4d3a-99da-71c5b1a22280
begin
	# ρ = minimum(ρvals)
	ρ = maximum(ρvals)
	mβ = mβ_vals[2]
	filter(f -> f[2]["ρ"] == ρ && f[2]["mβ"] == mβ, files)
	β2 = map(
		x -> x["β2"], 
		filter(f -> f[2]["ρ"] == ρ && f[2]["mβ"] == mβ, files) |> values
	) |> minimum
	
	params = getindex(
		files, 
		findfirst(f -> f["ρ"] == ρ && f["β2"] == β2 && f["mβ"] == mβ, files)
	)
end

# ╔═╡ 5efc4804-d5b9-4783-a6fe-9faf3dd89d40
f1 = readdlm(datdir * params["allele_freq_file"]);

# ╔═╡ 951c8684-cf03-4f0b-ab3a-46f3f361c56a
size(f1)

# ╔═╡ 6f77f00e-ddf2-4283-a0ee-6a72e1854173
L = Int(size(f1, 2)/2)

# ╔═╡ d76e1464-9f07-4c75-9ef6-e023c7f7882d
switch_info = let
	df = DataFrame(CSV.File(datdir * params["switch_info_file"]))
	df.time_to_next = vcat(df.t[2:end] - df.t[1:end-1], [missing])
	transform!(df, :h => (h -> map(x -> 1 - exp(-2*abs(x)/α), h)) => :β)
	transform!(df, 
		:h => 
		(h -> map(h) do h
			s = 2*abs(h)
			β = 1 - exp(-s/α)
			ν*β > x0 ? log(ν/(1-ν) * (β-x0)/x0)/s : missing
		end)
		=> :time_for_sweep
	)
	df
end

# ╔═╡ b708d3e2-7b36-4154-b5f8-632c9920d622
Δf_other_allele_sweep_time = map(eachrow(switch_info)) do r
	idx = Iterators.filter(!=(r.pos), 1:L-1)
	map(idx) do i
		if ismissing(r.time_for_sweep)
			missing
		elseif round(Int, r.t + r.time_for_sweep + 1) >= size(f1, 1)
			missing
		else
			f0 = f1[r.t + 1, 2*i]
			fend = f1[round(Int, r.t + r.time_for_sweep + 1), 2*i]
			Δf = if f0 < 0.05 || f0 > 0.95
				missing
			elseif (fend - f0) > 0
				(fend - f0) / r.β / (1-f0)
			elseif fend - f0 <= 0
				(fend - f0) / r.β / f0
			end
		end
	end
end

# ╔═╡ 04e37128-06dc-40f1-898a-633d185d2a5e
Δf_other_allele_time_to_next = map(eachrow(switch_info)) do r
	idx = Iterators.filter(!=(r.pos), 1:L-1)
	map(idx) do i
		if ismissing(r.time_to_next)
			missing
		elseif round(Int, r.t + r.time_to_next + 1) >= size(f1, 1)
			missing
		else
			f0 = f1[r.t + 1, 2*i]
			fend = f1[round(Int, r.t + r.time_to_next + 1), 2*i]
			Δf = if f0 < 0.05 || f0 > 0.95
				missing
			elseif (fend - f0) > 0
				(fend - f0) / r.β / (1-f0)
			elseif fend - f0 <= 0
				(fend - f0) / r.β / f0
			end
		end
	end
end

# ╔═╡ 19c25571-7736-4809-9192-7dda89b5f69f
dat_oa = let
	X = map(zip(
		switch_info.β, 
		Δf_other_allele_sweep_time, 
		Δf_other_allele_time_to_next
	)) do (β, df1, df2)
		idx = intersect(findall(!ismissing, df1), findall(!ismissing, df2))
		hcat(β * ones(length(idx)), df1[idx], df2[idx])
	end |> x -> vcat(x...)
	(
		β = Vector{Float64}(X[:,1]), 
		Δf_sweep_time = Vector{Float64}(X[:,2]), 
		Δf_time_to_next = Vector{Float64}(X[:,3]),
	)
end

# ╔═╡ cc663c4a-e860-4079-bac9-52e91036341d
let
	Δf_time_to_next = map(eachrow(switch_info)) do r
		i = r.h > 0 ? 2*r.pos-1 : 2*r.pos
		if ismissing(r.time_to_next)
			missing
		else
			t = round(Int, r.t + r.time_to_next)
			t < size(f1, 1) ? f1[t, i] : missing
		end
	end

	idx_fast = findall(eachrow(switch_info)) do r
		f = r.time_for_sweep < r.time_to_next
		ismissing(f) ? false : f
	end
	
	scatter(switch_info.β, Δf_time_to_next, label="")
	# scatter(switch_info.β[idx_fast], Δf_time_to_next[idx_fast], label="")
	plot!([0, 1], [0, 1], line=(:black), label="")

	plot!(
		xlabel = "β", 
		ylabel="Δf", 
		title = "Δf at start of next sweep\n <β>=$mβ, std(β)=$(round(sqrt(β2 - mβ^2), sigdigits=2)), ρ=$ρ",
		margin = 5mm,
	)
end

# ╔═╡ 10df71e1-9c67-4c4e-9cbd-b56e487aad3a
let
	Δf_time_for_sweep = map(eachrow(switch_info)) do r
		i = r.h > 0 ? 2*r.pos-1 : 2*r.pos
		if ismissing(r.time_for_sweep)
			missing
		else
			t = round(Int, r.t + r.time_for_sweep + 1)
			t < size(f1, 1) ? f1[t, i] - f1[r.t + 1, i] : missing
		end
	end

	scatter(switch_info.β, Δf_time_for_sweep, label="")
	plot!([1e-2, 1], [1e-2, ν], line=(:black), label="")

	plot!(
		xlabel = "β", 
		ylabel="Δf", 
		title = "Δf after 90%-time\n <β>=$mβ, std(β)=$(round(sqrt(β2 - mβ^2), sigdigits=2)), ρ=$ρ",
		margin = 5mm,
		# scale=:log10,
		# xlim=(1e-2, 1.5),
		# ylim=(1e-2, 1.5)
	)
end

# ╔═╡ c5456e2e-e725-4d07-904a-d0a94d5a69ac
let
	scatter(
		switch_info.time_to_next, switch_info.time_for_sweep;
		label="", marker=(2, 0.5), markerstrokewidth=0
	)
	plot!([1, 100], [1, 100], line=(:black, :dash), label="")
	plot!(
		xlabel = "Time to next sweep", 
		ylabel = "90%-time for sweep",
		title = "<β>=$mβ, std(β)=$(round(sqrt(β2 - mβ^2), sigdigits=2)), ρ=$ρ"
	)
end

# ╔═╡ 09e32347-9a28-49d2-96a2-ae57a647b2f0
let
	p = scatter(
		dat_oa.β, abs.(dat_oa.:Δf_sweep_time);
		color = 1, label="", markerstrokewidth=0, marker=(2.5, 0.5)
	)

	hline!([1], line=(:black, :dash), label="")
	# plot!(0.2:0.01:1, 1 ./(0.2:0.01:1), line=(:black, :dashdot), label="")

	# KDE smoothing
	β_smooth, Δf_oa_smooth = let
		βv = 0:0.01:1
		βv, smooth(dat_oa.β, abs.(dat_oa.:Δf_sweep_time), 0.05, βv)
	end
	plot!(β_smooth, Δf_oa_smooth, label="", color=:blue, line=(4, ))

	# 
	plot!(
		xlabel = "β", 
		ylabel="Δf (scaled)", 
		title = "Scaled Δf after 90%-time (other alleles)\n <β>=$mβ, std(β)=$(round(sqrt(β2 - mβ^2), sigdigits=2)), ρ=$ρ",
		margin = 5mm,
	)
end

# ╔═╡ 74434f0a-9a1f-454e-848b-fb3225a2d650
let
	p = scatter(
		dat_oa.β, abs.(dat_oa.Δf_time_to_next);
		color = 1, label="", markerstrokewidth=0, marker=(2.5, 0.5)
	)

	hline!([1], line=(:black, :dash), label="")
	# plot!(0.2:0.01:1, 1 ./(0.2:0.01:1), line=(:black, :dashdot), label="")

	# KDE smoothing
	β_smooth, Δf_oa_smooth = let
		βv = 0:0.01:1
		βv, smooth(dat_oa.β, abs.(dat_oa.Δf_time_to_next), 0.05, βv)
	end
	plot!(β_smooth, Δf_oa_smooth, label="", color=:blue, line=(4, ))

	# 
	plot!(
		xlabel = "β", 
		ylabel="Δf (scaled)", 
		title = "Scaled Δf at start of next sweep (other alleles)\n <β>=$mβ, std(β)=$(round(sqrt(β2 - mβ^2), sigdigits=2)), ρ=$ρ",
		margin = 5mm,
	)
end

# ╔═╡ 4d5f4de9-a442-47e9-a2c7-210610cf2fb8
md"""
Reason for the $1/\beta$ envelope curve: if an allele initially at $x_0$ has a positive frequency change after a partial sweep at another position, we expect $\Delta f = \beta(1-x_0)$. If the allele ends up fixating in the population due to other sweeps, we then observe $\Delta f^{max} = 1-x_0$. The ratio between these two quantities if $1/\beta$.
"""

# ╔═╡ a66bc3ad-6905-42d3-ad2b-0cddf40af9e7
md"## Tests"

# ╔═╡ 423ff8c5-c9d3-4a81-93ca-cba4188b01f4
md"""
For the first sweep, compare observed frequency to the analytical approximation, and whether the computed time is correct
"""

# ╔═╡ c5eed6d5-d5e5-49e7-aa2d-18bbc03783d3
let
	i = 20
	tvals = switch_info.t[i] .+ (0:300)
	pos = switch_info.pos[i]
	β = switch_info.β[i]
	s = 2*abs(switch_info.h[i])
	
	plot(tvals, f1[tvals .+ 1, 2*pos])
	hline!([ν * β])
	vline!([switch_info.time_for_sweep[i] + switch_info.t[i]])
	plot!(legend = :bottomright)
	
	plot!(tvals, map(t -> β*x0*exp(s*t) / (x0*exp(s*t) + β - x0), (tvals .- switch_info.t[i])))
end

# ╔═╡ Cell order:
# ╠═c7e4d226-f881-11ed-0717-0f30a6dc429a
# ╠═385ae3ab-f8d3-4d8d-9cca-43730702083c
# ╠═6a7e0c08-f9e0-475b-80f7-ebc2a69294a6
# ╠═167cdfe9-fa90-4982-92f0-430906199130
# ╠═b2ccf546-18aa-49a3-8499-ac21d4cc54ce
# ╠═37d544e6-e6a5-4e08-852c-d67820774fe5
# ╠═5efc4804-d5b9-4783-a6fe-9faf3dd89d40
# ╟─13826c8a-1770-4f6c-a5a3-6b0163fafa82
# ╠═fa04fcc4-081f-467f-b2ac-b890be897c83
# ╟─b8dea52e-be5c-47a3-9bc4-203fa360a035
# ╟─d76e1464-9f07-4c75-9ef6-e023c7f7882d
# ╟─60700a60-d523-4429-90f5-62ea201c3dfd
# ╟─2f25a0eb-96a7-4e87-acd0-cdc4c10e1b05
# ╟─cc663c4a-e860-4079-bac9-52e91036341d
# ╟─10df71e1-9c67-4c4e-9cbd-b56e487aad3a
# ╟─c5456e2e-e725-4d07-904a-d0a94d5a69ac
# ╟─8a929898-b50e-43d3-bf14-bf8e72ca8ee2
# ╠═951c8684-cf03-4f0b-ab3a-46f3f361c56a
# ╟─b58c5cab-b329-4207-9f05-4e4ce0e9e98f
# ╠═6f77f00e-ddf2-4283-a0ee-6a72e1854173
# ╟─b708d3e2-7b36-4154-b5f8-632c9920d622
# ╟─04e37128-06dc-40f1-898a-633d185d2a5e
# ╠═19c25571-7736-4809-9192-7dda89b5f69f
# ╠═1fe37cb1-b70a-4d3a-99da-71c5b1a22280
# ╟─09e32347-9a28-49d2-96a2-ae57a647b2f0
# ╟─74434f0a-9a1f-454e-848b-fb3225a2d650
# ╟─4d5f4de9-a442-47e9-a2c7-210610cf2fb8
# ╟─a66bc3ad-6905-42d3-ad2b-0cddf40af9e7
# ╟─423ff8c5-c9d3-4a81-93ca-cba4188b01f4
# ╠═c5eed6d5-d5e5-49e7-aa2d-18bbc03783d3
