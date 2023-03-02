### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 2fe14f34-b81c-11ed-3729-67545b356a8b
begin
	using Pkg; Pkg.activate("../../")
	using Chain
	using EasyFFTs
	using LinearAlgebra
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR # prefix functions with PSS, e.g. `PSS.Parameters`
	using StatsBase
end

# ╔═╡ 2754bfbd-9e4c-4881-a7a9-c5f4257141d1
PlutoUI.TableOfContents()

# ╔═╡ 1047ab91-66b9-4e95-8b1c-c98e5ee6d48d
md"# Parameters"

# ╔═╡ 4e258ed8-a888-407f-ad22-c4362495e03a
begin
	α = 3
	δ = 1
end;

# ╔═╡ fb9110ed-b998-4d38-afe0-230c60b37872
md"""
Setting up cross-immunity: let us keep the equilibrium frequency $x$ constant and vary the distance $d$ of cross-immunity terms from $1$: 

$$\begin{aligned}
x &= \frac{1-f}{(1-f) + (1-b)} \;\; \text{(equilibrium frequency)}\\
& \\
d &= \frac{(1-b) + (1-f)}{2} \;\; \text{(distance from $1$ for $b$ and $f$)}\\ 
\end{aligned}$$

Inverting this gives 

$$f = 1-2xd, \; f = 1-2(1-x)d.$$
"""

# ╔═╡ 59a1579d-d0c1-4378-bb43-f85b568483d8
x_slider = @bind x Slider(0.05:0.05:0.95, show_value=true, default = 0.65)

# ╔═╡ d60e523f-be75-48d5-9f00-6ea0681f0a33
d_slider = @bind d Slider(0.025:0.025:0.95, show_value=true, default = 0.2)

# ╔═╡ 1d63d866-7ec0-4c59-bc50-a301a4324384
K = begin
	f = 1-2*x*d
	b = 1 - 2*(1-x)*d
	[1 b; f 1]
end

# ╔═╡ 4e012216-6cea-471f-ba57-ba7daeb77286
γ_slider = let
	γmax = 4*(α-δ)/α^2*(1-b)*(1-f)/(1-b*f)
	r = 10 .^ range(-4, log10(γmax - 1e-6), length=25)
	_γ = @bind γ Slider(r, show_value=true, default = 0.005)
	md"γ = $(_γ)"
end

# ╔═╡ a451e995-e1e4-4845-a180-0072faf864ae
params = let
	N = 2
	M = 1
	PSS.Parameters(; N, M, α, γ, δ)
end

# ╔═╡ 24cdb593-8465-4180-a42a-082bcb6a01a8
# Simulation time
T = 5/γ

# ╔═╡ 0fb3d755-b409-47d0-bc4d-75a002ba96d6
md"""
**Note**: One can write `x_slider` or `d_slider` in any cell to tune the sliders
"""

# ╔═╡ b2148183-228f-45bf-a879-b5e6039ea5f4
md"# Simulation"

# ╔═╡ 7f8d658c-c329-4092-87ae-cfe33d44237e
md"Defining a region, initially with no mutant"

# ╔═╡ 599065eb-0f6d-4d70-85ff-8e1ee6e2ad3f
region = let
	S0 = .4
	I0 = 1e-6
	C0 = 1e-9
	R0 = 1 - S0 - I0 - C0
	wt = PSS.Virus(; S=S0, I=I0, C=C0, R=R0)
	mut = PSS.Virus(; S=1, I=0, C=0., R=0.)
	PSS.Region(; viruses=[wt, mut], K)
end

# ╔═╡ eb246d46-a899-4330-b430-a1c0ff5d8bc7
md"""
Simulating for time $T$ without the mutant to reach the initial equilibrium. Then, introduce a small quantity of the mutant in the region to obtain `state_init`. 
"""

# ╔═╡ 21040858-5c97-43a3-b38f-4e4ce5c90ef7
state_init = let
	state = PSS.SIRState(; regions=[region], parameters=params)
	sol = PSS.simulate(state, (0, T))
	PSS.set_infected(sol(T), 2, 1e-6) # infect with mutant
end;

# ╔═╡ 29993cd9-f004-49fa-88da-82b436373c6a
md"Finally, simulate again for time $T$."

# ╔═╡ 8abbc0cd-e809-4d70-81d6-e4f783c94b53
sol = PSS.simulate(state_init, (0, T)); # simulated solution

# ╔═╡ 2e9996f8-093b-4d9d-b00b-b1d85c2a3267
md"Extract the final state, and also compute the analytical equilibrium."

# ╔═╡ 5398dc6c-5847-43a7-b0ef-2d4f7353a1f3
state_final = sol(T);

# ╔═╡ 05f3e103-28a2-466e-829f-ea302f921712
state_eq = PSS.equilibrium(state_init); # analytical equilibrium

# ╔═╡ 3183ebc7-d067-4b8e-b7a7-1ff619b0f7d7
md"# Oscillations"

# ╔═╡ ac9a319b-f7e4-418d-aefa-9de2e3fd53dd
md"## Theory"

# ╔═╡ 806faa77-c2ee-4af6-8dcd-e8009a5dbf26
md"""
Around equilibrium, the solution behaves as 

$$\dot{X} = QX,$$

where $X = [S^{wt}, S^m, I^{wt}, I^m]$ and 

$$\begin{aligned}
	Q &= \begin{pmatrix}
	-\alpha\gamma/\delta & 0 & -\delta & -\delta b\\
	0 & -\alpha\gamma/\delta & -\delta f & -\delta\\
	g_1 & 0 & 0 & 0\\
	0 & g_2 & 0 & 0\\
	\end{pmatrix}\\
	&\\
	g_1 & = \frac{\gamma}{\delta}(\alpha-\delta)\frac{1-b}{1-bf}\\
	g_2 & = \frac{\gamma}{\delta}(\alpha-\delta)\frac{1-f}{1-bf}.\\
\end{aligned}$$
"""

# ╔═╡ a4131661-44ac-48c4-8044-017dbcbf1c42
md"""
For low enough $\gamma$, we find the eigenvalues to be 

$$\begin{aligned}
	\lambda_1 &= -\frac{1}{2}\frac{\alpha}{\delta}\gamma \pm i\left( \gamma(\alpha-\delta)  - \frac{1}{4}\frac{\alpha^2\gamma^2}{\delta^2}\right)^{1/2}\\
	&\\
	\lambda_2 &= -\frac{1}{2}\frac{\alpha}{\delta}\gamma \pm i\left( \gamma(\alpha-\delta)\frac{(1-b)(1-f)}{1-bf}  - \frac{1}{4}\frac{\alpha^2\gamma^2}{\delta^2}\right)^{1/2}.\\
\end{aligned}$$

*Note*: the above decomposition is only valid if the expressions inside the squared roots are positive. We thus have the following condition on $\gamma$ for this to be valid:

$$\gamma \leq \frac{4(\alpha-\delta)\delta^2}{\alpha^2}\frac{(1-b)(1-f)}{1-bf}.$$

This will always be true in our setting since we assume $\gamma \ll \alpha, \delta$ and we do not have $b,f\ll 1$.

"""

# ╔═╡ b57d0235-2e8d-4bda-b080-0bff5b993104
md"""
We thus have the following behavior close to equilibrium
- the rate at which the return to equilibrium happens after a perturbation is $\frac{\alpha\gamma}{2\delta}$. Note that the equilibrium is stable;
- there are two oscillatory frequencies, simplified in the small $\gamma$ limit: 
$$\begin{aligned}
	\nu_{high} &= \frac{\sqrt{\gamma(\alpha-\delta)}}{2\pi},\\
	\nu_{low} &= \nu_{high}\sqrt{\frac{(1-b)(1-f)}{1-bf}}.\\
\end{aligned}$$
"""

# ╔═╡ ee492499-0346-4f97-a44d-fb0196e6f6cf
let
	lw = 2
	
	γmax = 4*(α-δ)/α^2*(1-b)*(1-f)/(1-b*f)
	γvals = range(1e-5, γmax - 1e-6, length=1000)

	ν_high(γ) = sqrt(γ*(α-δ) - (α*γ/2/δ)^2) / 2/3.14
	νs_high(γ) = sqrt(γ*(α-δ)) / 2/3.14

	ν_low(γ) = sqrt(γ*(α-δ)*(1-b)*(1-f)/(1-b*f) - (α*γ/2/δ)^2) / 2/3.14
	νs_low(γ) = νs_high(γ)*sqrt((1-b)*(1-f)/(1-b*f))

	p = plot(
		xlabel = "γ",
		title = "Oscillation frequencies",
		xscale = :log10,
		yscale = :log10
	)
	# high
	plot!(γvals, νs_high.(γvals), label="ν high")
	plot!(
		γvals, ν_high.(γvals); 
		line=(:dash), color=1, label = "ν high (exact)"
	)

	# low
	plot!(γvals, νs_low.(γvals), label="ν low", color=2)
	plot!(
		γvals, ν_low.(γvals); 
		line=(:dash), color=2, label = "ν low (exact)"
	)

	# upper threshold for gamma
	vline!([γmax], label="", line=(:black))
end

# ╔═╡ cb7b124d-e6ff-482a-af8b-054502fefad7
begin
	λ1 = -α/δ*γ/2 + im * sqrt(γ*(α-δ) - (α/δ*γ)^2/4)
	λ2 = -α/δ*γ/2 + im * sqrt(γ*(α-δ)*(1-b)*(1-f)/(1-b*f) - (α/δ*γ)^2/4)
end;

# ╔═╡ d83ea198-be76-4b6f-af80-87f1ae8a1fae
begin
	ν_high = sqrt(γ*(α-δ)) / 2/3.14
	ν_low = ν_high * sqrt((1-b)*(1-f)/(1-b*f))
end;

# ╔═╡ 0f976fde-3829-4aa4-8047-f5645efa6d03
md"## Comparison with simulations"

# ╔═╡ bae05ba2-fcc4-4a53-8a27-0b5603c823fc
md"""
To compare with simulations, we will compute the FFT of the time series of different quantities ($I$, $S$, frequency of mutant...) in a time interval $[2\gamma^{-1}, T]$ where $T$ is large w.r. to $\gamma^{-1}$. Essentially, we expect our theory to be valid around equilibrium. 

*Note*: If one computes the FFT in a time interval $[0, T]$, it does not work as well. 
"""

# ╔═╡ 2344d4b7-e50d-4bcc-adb7-9f8f3f58b88b
md"""
*Note*: For the FFT to look a bit nicer, we work with scaled quantities. For example, considering the time series $I^{wt}(t)$, we rescale it to 

$$I^{wt}_{scaled}(t) = I^{wt}(t) - \frac{1}{\Delta t}\int_{t_0}^{t_0+\Delta t} I^{wt}(t)\text{d}t$$
"""

# ╔═╡ e83a818b-f738-452c-b24e-d544a1f904c3
fs = 1 # sampling frequency for fft

# ╔═╡ 9d81f6f3-b341-4a56-ac8c-c63de91c541f
tvals = range(sol.tspan..., step=1/fs) # relevant range of time values 

# ╔═╡ 9b46ca07-f06d-4fe3-aecf-4fa0fee14939
tmin_fft = 2/γ

# ╔═╡ 53f51316-f10d-4cbb-a9b7-b6f392acb4d0
tvals_fft = range(tmin_fft, T, step = 1/fs)

# ╔═╡ cb23bd65-1b98-4242-9a08-ae2c4fb4c208
md"### Infectious"

# ╔═╡ 050fabd1-7045-42ba-a979-e7219d61dc98
Ivals_scaled = let
	Iwt = sol[tvals_fft, 1, 1, :I]
	Im = sol[tvals_fft, 1, 2, :I]
	Itot = Iwt + Im

	(
		Iwt = Iwt .- mean(Iwt),
		Im = Im .- mean(Im),
		Itot = Itot .- mean(Itot)
	)
end;

# ╔═╡ 383a27dd-c361-4d91-a69c-c90dc9f5ee61
let
	p1 = plot(
		xlabel="time", 
		ylabel="",
		title="Infectious (scaled)",
		legend=:bottomright,
	)
	

	plot!(tvals_fft, Ivals_scaled.Iwt, label="wt")
	plot!(tvals_fft, Ivals_scaled.Im, label="mutant")

	p2 = plot(
		xlabel="time", 
		ylabel="",
		title="Infectious (scaled)",
		legend=:bottomright,
	)
	plot!(tvals_fft, Ivals_scaled.Itot, label="total", line=(:black))
	
	plot(p1, p2, layout=grid(1,2), size = (900, 450))
end

# ╔═╡ 2fb223fe-f2ef-4f47-9713-d37c8da4261e
let
	ps = map((:Iwt, :Im, :Itot)) do a
		I = getproperty(Ivals_scaled, a)
	
		fft = easyfft(I, fs)
		
		p = plot(
			title = "FFT $a",
			xlim = (-0.01, 5*ν_high),
			xlabel = "frequency"
		)
		plot!(fft.freq, magnitude(fft), label="")
		vline!([ν_low], label="ν low")
		vline!([ν_high], label="ν high")
	end
	plot(ps..., layout = grid(1,3), size = (900,300))
end

# ╔═╡ a0089de5-2de2-4b0f-a57c-29e177d8696c
md"### Susceptibles"

# ╔═╡ 244bed9f-9bff-496e-94a0-2194fa94a2a5
Svals_scaled = let
	Swt = sol[tvals_fft, 1, 1, :S]
	Sm = sol[tvals_fft, 1, 2, :S]
	Stot = Swt + Sm

	(
		Swt = Swt .- mean(Swt),
		Sm = Sm .- mean(Sm),
		Stot = Stot .- mean(Stot)
	)
end;

# ╔═╡ 30a996df-1ed0-47e7-be32-c710fdaecf7c
let
	p1 = plot(
		xlabel="time", 
		ylabel="",
		title="Susceptible (scaled)",
		legend=:bottomright,
	)
	

	plot!(tvals_fft, Svals_scaled.Swt, label="wt")
	plot!(tvals_fft, Svals_scaled.Sm, label="mutant")

	p2 = plot(
		xlabel="time", 
		ylabel="",
		title="Susceptible (scaled)",
		legend=:bottomright,
	)
	plot!(tvals_fft, Svals_scaled.Stot, label="total", line=(:black))
	
	plot(p1, p2, layout=grid(1,2), size = (900, 450))
end

# ╔═╡ 7b54e7c1-f4fc-4689-a96d-66ff411f279d
let
	ps = map((:Swt, :Sm, :Stot)) do a
		S = getproperty(Svals_scaled, a)
	
		fft = easyfft(S, fs)
		
		p = plot(
			title = "FFT $a",
			xlim = (-0.01, 5*ν_high),
			xlabel = "frequency"
		)
		plot!(fft.freq, magnitude(fft), label="")
		vline!([ν_low], label="ν low")
		vline!([ν_high], label="ν high")
	end
	plot(ps..., layout = grid(1,3), size = (900,300))
end

# ╔═╡ 0c7373d7-db48-4855-b164-67a4757076f4
md"### Frequencies"

# ╔═╡ 4cd4f525-adff-402a-b51a-3b02c3afd1a2
freq_scaled = let
	Iwt = sol[tvals_fft, 1, 1, :I]
	Im = sol[tvals_fft, 1, 2, :I]
	freq = Im ./ (Im + Iwt)
	freq .- mean(freq)
end;

# ╔═╡ c3c3f7ce-3083-4717-860d-7c5afbc12841
let
	p = plot(
		xlabel="time", 
		ylabel="",
		title="Frequency of mutant (scaled)",
		legend=:bottomright,
	)
	plot!(tvals_fft, freq_scaled, label="")
end

# ╔═╡ 98f9a526-b28e-440e-b772-131a21cae485
let
	fft = easyfft(freq_scaled, fs)
	
	p = plot(
		title = "FFT of mutant frequency",
		xlim = (-0.01, 5*ν_high),
		xlabel = "frequency"
	)
	plot!(fft.freq, magnitude(fft), label="")
	vline!([ν_low], label="ν low")
	vline!([ν_high], label="ν high")
end

# ╔═╡ Cell order:
# ╠═2fe14f34-b81c-11ed-3729-67545b356a8b
# ╠═2754bfbd-9e4c-4881-a7a9-c5f4257141d1
# ╟─1047ab91-66b9-4e95-8b1c-c98e5ee6d48d
# ╠═4e258ed8-a888-407f-ad22-c4362495e03a
# ╟─4e012216-6cea-471f-ba57-ba7daeb77286
# ╠═a451e995-e1e4-4845-a180-0072faf864ae
# ╠═24cdb593-8465-4180-a42a-082bcb6a01a8
# ╟─fb9110ed-b998-4d38-afe0-230c60b37872
# ╠═59a1579d-d0c1-4378-bb43-f85b568483d8
# ╠═d60e523f-be75-48d5-9f00-6ea0681f0a33
# ╠═1d63d866-7ec0-4c59-bc50-a301a4324384
# ╟─0fb3d755-b409-47d0-bc4d-75a002ba96d6
# ╟─b2148183-228f-45bf-a879-b5e6039ea5f4
# ╟─7f8d658c-c329-4092-87ae-cfe33d44237e
# ╠═599065eb-0f6d-4d70-85ff-8e1ee6e2ad3f
# ╟─eb246d46-a899-4330-b430-a1c0ff5d8bc7
# ╠═21040858-5c97-43a3-b38f-4e4ce5c90ef7
# ╟─29993cd9-f004-49fa-88da-82b436373c6a
# ╠═8abbc0cd-e809-4d70-81d6-e4f783c94b53
# ╟─2e9996f8-093b-4d9d-b00b-b1d85c2a3267
# ╠═5398dc6c-5847-43a7-b0ef-2d4f7353a1f3
# ╠═05f3e103-28a2-466e-829f-ea302f921712
# ╟─3183ebc7-d067-4b8e-b7a7-1ff619b0f7d7
# ╟─ac9a319b-f7e4-418d-aefa-9de2e3fd53dd
# ╟─806faa77-c2ee-4af6-8dcd-e8009a5dbf26
# ╟─a4131661-44ac-48c4-8044-017dbcbf1c42
# ╟─b57d0235-2e8d-4bda-b080-0bff5b993104
# ╟─ee492499-0346-4f97-a44d-fb0196e6f6cf
# ╠═cb7b124d-e6ff-482a-af8b-054502fefad7
# ╠═d83ea198-be76-4b6f-af80-87f1ae8a1fae
# ╟─0f976fde-3829-4aa4-8047-f5645efa6d03
# ╟─bae05ba2-fcc4-4a53-8a27-0b5603c823fc
# ╟─2344d4b7-e50d-4bcc-adb7-9f8f3f58b88b
# ╠═e83a818b-f738-452c-b24e-d544a1f904c3
# ╠═9d81f6f3-b341-4a56-ac8c-c63de91c541f
# ╠═9b46ca07-f06d-4fe3-aecf-4fa0fee14939
# ╠═53f51316-f10d-4cbb-a9b7-b6f392acb4d0
# ╟─cb23bd65-1b98-4242-9a08-ae2c4fb4c208
# ╠═050fabd1-7045-42ba-a979-e7219d61dc98
# ╟─383a27dd-c361-4d91-a69c-c90dc9f5ee61
# ╟─2fb223fe-f2ef-4f47-9713-d37c8da4261e
# ╟─a0089de5-2de2-4b0f-a57c-29e177d8696c
# ╠═244bed9f-9bff-496e-94a0-2194fa94a2a5
# ╟─30a996df-1ed0-47e7-be32-c710fdaecf7c
# ╟─7b54e7c1-f4fc-4689-a96d-66ff411f279d
# ╟─0c7373d7-db48-4855-b164-67a4757076f4
# ╠═4cd4f525-adff-402a-b51a-3b02c3afd1a2
# ╟─c3c3f7ce-3083-4717-860d-7c5afbc12841
# ╟─98f9a526-b28e-440e-b772-131a21cae485
