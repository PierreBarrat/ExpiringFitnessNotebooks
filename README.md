# ExpiringFitnessNotebooks

## Install

- install [julia v1.8.5](https://julialang.org/downloads/). *Note*: [juliaup](https://github.com/JuliaLang/juliaup) makes managing julia versions easier.
- `git clone` this somewhere. Navigate to the repository. 
- start julia and enter the following

```
using Pkg; 
Pkg.activate(".") # activates the local environment
Pkg.resolve() # create the Manifest.toml file
Pkg.instantiate() # download and precompile all packages
```

  Note that precompilation could take some time, especially `PartialSweepSIR` which depends on `OrdinaryDiffEq`.

- start `Pluto` by typing `using Pluto; Pluto.run()`. Once the Pluto server is ready, use the explorer in the Pluto window to open any notebook. The first run of each notebook can be quite long. 

*Note*: Some hints about how to use Pluto can be found on the [github page](https://github.com/fonsp/Pluto.jl), in the [wiki/faq](https://github.com/fonsp/Pluto.jl/wiki) and in the featured notebooks shown on the welcome screen. 




