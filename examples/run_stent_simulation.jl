# =========================================================================
# examples/run_stent_simulation.jl
# Demonstration of a 30-day degradation profile for a biodegradable cardiovascular stent.
# =========================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, "..")) # Automatically activates the local environment

using CorrosionPredictor
using Unitful

println("==================================================")
println("🚀 Initializing Cardiovascular Stent Simulation...")
println("==================================================")

# 1. Define the physiological environment for the stent
# Including all necessary transport and boundary layer parameters
env_params = Dict(
    "pH_local" => 7.4,
    "T_temperature" => 310.15u"K",
    "c_O2_bulk" => 0.2u"mol/m^3",
    "c_Zn2_bulk" => 0.0u"mol/m^3",
    "c_H2O2_bulk" => 0.0u"mol/m^3",
    "D_O2" => 2.0e-9u"m^2/s",
    "D_Zn2" => 0.7e-9u"m^2/s",
    "D_H2O2" => 1.4e-9u"m^2/s",
    "boundary_layer_thickness" => 10.0e-6u"m",
    "initial_mass" => 15.0,
    "dt_initial" => 3600.0,
    "mechanical_rupture_days" => 5.0,
    
    # --- Kinetic Parameters (Triggers Dynamic Butler-Volmer Solver) ---
    "i0_zinc" => 1.0e-2u"A/m^2",    # Exchange current density for Zinc
    "alpha_a" => 0.5,               # Anodic charge transfer coefficient
    "alpha_c" => 0.5,               # Cathodic charge transfer coefficient
    "ionic_strength_M" => 0.15      # Physiological saline
)

# 2. Run the decoupled spatiotemporal simulation for 30 days
println("Material: ", ZINC.name)
println("Simulated duration: 30 days")
println("Mode: Dynamic Mixed-Potential Kinetics Active")

t_final = 86400.0 * 30 
history = run_corrosion_simulation(ZINC, env_params, t_final)

# 3. Output results
println("\n✅ Simulation complete!")
println("Final Mass: ", round(history.mass[end], digits=2), " mg")
println("Final Oxide Layer Thickness: ", round(history.layer[end], digits=2), " μm")
println("Final Interface pH: ", round(history.pH[end], digits=2))