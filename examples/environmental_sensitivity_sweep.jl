# =========================================================================
# examples/environmental_sensitivity_sweep.jl
# Demonstration of a physiological parameter sweep.
# Evaluates how different bulk oxygen concentrations (hypoxia vs normoxia)
# impact the dynamic mixed potential and the corrosion rate.
# =========================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CorrosionPredictor
using Unitful
using Plots

function run_sweep()
    println("==================================================")
    println("🔬 Running Environmental Sensitivity Sweep...")
    println("==================================================")

    # Define baseline physiological parameters
    base_env = Dict(
        "pH_local" => 7.4,
        "T_temperature" => 310.15u"K",
        "c_Zn2_bulk" => 0.0u"mol/m^3",
        "c_H2O2_bulk" => 0.0u"mol/m^3",
        "D_O2" => 2.0e-9u"m^2/s",
        "D_Zn2" => 0.7e-9u"m^2/s",
        "D_H2O2" => 1.4e-9u"m^2/s",
        "boundary_layer_thickness" => 10.0e-6u"m",
        "initial_mass" => 15.0,
        "dt_initial" => 3600.0,
        "mechanical_rupture_days" => 0.0,
        "i0_zinc" => 1.0e-2u"A/m^2",
        "alpha_a" => 0.5,
        "alpha_c" => 0.5,
        "ionic_strength_M" => 0.15
    )

    # Define a sweep of oxygen concentrations (from severe hypoxia to normoxia)
    o2_levels = [0.02, 0.05, 0.10, 0.15, 0.22] .* u"mol/m^3"
    final_masses = Float64[]
    final_layer_thicknesses = Float64[]

    t_final = 86400.0 * 30 # Simulate 30 days for each scenario

    for c_O2 in o2_levels
        # Update environment for current iteration
        current_env = copy(base_env)
        current_env["c_O2_bulk"] = c_O2
        
        # Run the decoupled simulation
        history = run_corrosion_simulation(ZINC, current_env, t_final)
        
        # Extract end-state metrics
        push!(final_masses, history.mass[end])
        push!(final_layer_thicknesses, history.layer[end])
        
        println("Tested [O2] = ", round(typeof(1.0u"mol/m^3"), c_O2, digits=2), 
                " | Remaining Mass: ", round(history.mass[end], digits=2), " mg")
    end

    # Create results directory if it does not exist
    results_dir = joinpath(@__DIR__, "..", "results")
    mkpath(results_dir)

    # Generate a visual artifact of the sweep
    p = plot(ustrip.(o2_levels), final_masses, 
             title="Stent Mass after 30 Days vs Bulk Oxygen",
             xlabel="Bulk O2 Concentration (mol/m³)", 
             ylabel="Remaining Mass (mg)",
             marker=:circle, lw=2, legend=false, color=:blue, grid=true)

    # Save the plot in the new results folder
    output_path = joinpath(results_dir, "environmental_sweep_results.png")
    savefig(p, output_path)

    println("\n✅ Sweep complete! Artifact saved as '", output_path, "'")
end

# Execute the simulation
run_sweep()