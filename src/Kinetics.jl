# =========================================================================
# src/Kinetics.jl
# Master Router and Explicit Submodule for the Electrochemical Kinetics Domain.
# Encapsulated in an isolated namespace to prevent global scope contamination.
# =========================================================================

module Kinetics

using Unitful
using Roots # Explicit dependency for root finding algorithms

# Import core material, dimensional, and transport/physics dependencies cleanly from the parent module context
import ..CorrosionPredictor: BioresorbableMaterial, Material, ZINC, MAGNESIUM,
                             CorrosionRate, Concentration, Length, Area, MolarFlux,
                             Voltage, CurrentDensity, GasConstant, FaradayConstant, Temperature,
                             calculate_interface_concentration, calculate_film_resistance_factor

# Export the public electrochemical API of this subdomain
export calculate_activity_coefficient, calculate_nernst_potential
export calculate_butler_volmer_flux, corrosion_rate_to_flux
export find_mixed_potential, calculate_stoichiometric_fluxes, calculate_surface_environment

# Splitting the module into cohesive, isolated sub-files within the Kinetics namespace
include("Kinetics/Thermodynamics.jl")
include("Kinetics/FaradaicFlux.jl")
include("Kinetics/MixedPotentialSolver.jl")

end # module Kinetics