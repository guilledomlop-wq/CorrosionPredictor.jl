# CorrosionPredictor.jl 🧬⚙️
Predictive corrosion simulation platform built in Julia.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/guilledomlop-wq/CorrosionPredictor.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/guilledomlop-wq/CorrosionPredictor.jl/actions)

**CorrosionPredictor.jl** is a high-performance, spatiotemporal simulation engine designed for the in-silico study of advanced cardiovascular biodegradable stents. This platform bridges rigorous technical simulation stages with high-level scientific milestones. 

Written purely in Julia, it leverages multiple dispatch, strictly validated dimensional types (`Unitful.jl`), and symbolic PDE modeling (`ModelingToolkit.jl`) to decouple fast electrochemical kinetics from slow macroscopic transport.

## 🚀 Key Architectural Features

1. **Dimensional Type Safety (`Types.jl` & `Materials.jl`):** A centralized type system ensures that all SciML solvers and physical equations are protected against dimensional mismatches. Materials are implemented using abstract types and multiple dispatch for zero-overhead extensibility.
2. **Split-Operator Kinetics (`Kinetics.jl`):** An isolated submodule that calculates mixed-potential electrochemical roots, handling Butler-Volmer fluxes, dynamic Nernst potential shifts, and dynamic diffusion limits without polluting the global namespace.
3. **Symbolic Transport (`Transport.jl`):** Generates 1D spatio-temporal diffusion PDE systems using the Method of Lines (`MethodOfLines.jl`), ensuring mathematically robust mass-transport simulations.
4. **Adaptive Simulation Engine (`Simulator.jl`):** Features adaptive time-stepping with asymptotic recovery and mechanical rupture checks to simulate the physical pulsing of arteries on the stent's brittle protective layer.

## 📦 Installation

This package is currently in a pre-release state for portfolio showcase. You can add it locally or clone it:

```julia
julia> ]
pkg> add [https://github.com/guilledomlop-wq/CorrosionPredictor.jl](https://github.com/guilledomlop-wq/CorrosionPredictor.jl)
🔬 Quick Start
Julia
using CorrosionPredictor
using Unitful

# 1. Load biological environment parameters
env_params = Dict(
    "pH_local" => 7.4,
    "T_temperature" => 310.15u"K",
    "c_O2_bulk" => 0.2u"mol/m^3",
    "initial_mass" => 15.0 # mg
)

# 2. Run the decoupled spatiotemporal simulation
history = run_corrosion_simulation(ZINC, env_params, 86400.0 * 30) # 30 days
🛠️ Testing
The package includes a comprehensive test suite validating the physical bounds, unit consistency, and thermodynamic algorithms.

Julia
pkg> test