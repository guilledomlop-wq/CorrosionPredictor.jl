# =========================================================================
# src/CorrosionPredictor.jl
# Main entry point and orchestrator for the CorrosionPredictor module.
# =========================================================================

module CorrosionPredictor

using Unitful
using ModelingToolkit, MethodOfLines, DomainSets

# =========================================================================
# MODULE INCLUSION ORDER (Dependency Cascade)
# =========================================================================

include("Types.jl")
include("Materials.jl")
include("PhysicsCore.jl")

include("Kinetics.jl")
using .Kinetics      

include("Physics.jl")
include("Transport.jl")

# 6. Main Simulator
include("Simulator.jl")

# =========================================================================
# EXPORTS
# =========================================================================

# --- Types ---
export Length, Area, FluidVolume, Density, MolarMass, Concentration, DiffusionCoefficient, MolarFlux, CorrosionRate, MassRate, Voltage, CurrentDensity, GasConstant, FaradayConstant, Temperature, CapacitancePerArea
export BioresorbableMaterial, Material, ZincType, MagnesiumType, ZINC, MAGNESIUM

# --- Core Physics ---
export calculate_film_resistance_factor, calculate_interface_concentration
export calculate_mass_loss_rate, calculate_layer_growth_rate

# --- Kinetics ---
export calculate_activity_coefficient, calculate_nernst_potential, calculate_butler_volmer_flux, corrosion_rate_to_flux, find_mixed_potential, calculate_stoichiometric_fluxes, calculate_surface_environment

# --- Transport ---
export calculate_diffusive_flux, build_diffusion_model

# --- Main Simulator ---
export run_corrosion_simulation

end # module