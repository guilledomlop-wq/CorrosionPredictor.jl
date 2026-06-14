# =========================================================================
# examples/04_mechanical_rupture_effect.jl
# Demonstration of the mechanical rupture effect on the protective oxide layer
# caused by the physiological pulsing of the artery.
# Compares a static tissue environment vs a pulsing biomechanical environment.
# =========================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CorrosionPredictor
using Unitful

println("==================================================")
println("💓 Simulating Biomechanical Artery Pulsing vs Static...")
println("==================================================")

# Define the base physiological environment
base_env = Dict(
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
    "i0_zinc" => 1.0e-2u"A/m^2",
    "alpha_a" => 0.5,
    "alpha_c" => 0.5,
    "ionic_strength_M" => 0.15
)

t_final = 86400.0 * 30 # 30 days

# --- Scenario 1: Static Tissue (No Rupture) ---
# The oxide layer grows uninterrupted, maximizing physical protection
env_static = copy(base_env)
env_static["mechanical_rupture_days"] = 0.0
history_static = run_corrosion_simulation(ZINC, env_static, t_final)

# --- Scenario 2: Pulsing Artery (Rupture every 7 days) ---
# The biomechanical pulsing destroys 80% of the brittle protective layer periodically
env_pulsing = copy(base_env)
env_pulsing["mechanical_rupture_days"] = 7.0
history_pulsing = run_corrosion_simulation(ZINC, env_pulsing, t_final)

# Output comparative results
println("\n--- Scenario 1: Static Tissue (Perfect Protective Layer) ---")
println("Final Mass:                  ", round(history_static.mass[end], digits=2), " mg")
println("Final Oxide Layer Thickness: ", round(history_static.layer[end], digits=2), " μm")

println("\n--- Scenario 2: Pulsing Artery (Mechanical Rupture) ---")
println("Final Mass:                  ", round(history_pulsing.mass[end], digits=2), " mg")
println("Final Oxide Layer Thickness: ", round(history_pulsing.layer[end], digits=2), " μm")

println("\n✅ Biomechanical comparison complete!")