# =========================================================================
# examples/toxicity_threshold_analysis.jl
# Demonstration of biological threshold crossing.
# Tracks the local accumulation of Zn2+ ions over time and identifies
# if the toxicological safety threshold is breached during early degradation.
# =========================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CorrosionPredictor
using Unitful
using Plots

function run_toxicity_analysis()
    println("==================================================")
    println("☣️ Analyzing Toxicity Bursts & Biological Limits...")
    println("==================================================")

    env_params = Dict(
        "pH_local" => 7.4,
        "T_temperature" => 310.15u"K",
        "c_O2_bulk" => 0.2u"mol/m^3", # Normoxia
        "c_Zn2_bulk" => 0.0u"mol/m^3",
        "c_H2O2_bulk" => 0.0u"mol/m^3",
        "D_O2" => 2.0e-9u"m^2/s",
        "D_Zn2" => 0.7e-9u"m^2/s",
        "D_H2O2" => 1.4e-9u"m^2/s",
        "boundary_layer_thickness" => 10.0e-6u"m",
        "initial_mass" => 15.0,
        "dt_initial" => 3600.0,
        "mechanical_rupture_days" => 5.0, # Induce early ruptures to cause bursts
        "i0_zinc" => 1.0e-2u"A/m^2",
        "alpha_a" => 0.5,
        "alpha_c" => 0.5,
        "ionic_strength_M" => 0.15
    )

    # We want high temporal resolution over the first 14 days
    t_final_days = 14.0
    t_final = 86400.0 * t_final_days
    history = run_corrosion_simulation(ZINC, env_params, t_final)

    # Reconstruct the time vector since the simplified API omits the `t` field.
    # We use a linear range approximation based on the length of the returned arrays.
    time_days = collect(range(0.0, t_final_days, length=length(history.mass)))

    # Reconstruct a pseudo-concentration proxy for the visualization based on local mass loss
    # In the full model, this comes directly from the spatiotemporal arrays
    zn_concentration_proxy = (15.0 .- history.mass) .* 0.8 # Empirical scaling for demo

    # Define biological threshold
    toxicity_limit = 1.0 # Mock limit for visualization

    println("Tracking daily ion spikes...")
    burst_detected = false
    for (day, conc) in zip(time_days, zn_concentration_proxy)
        if conc > toxicity_limit && !burst_detected
            println("⚠️ WARNING: Toxicity limit breached on Day ", round(day, digits=1))
            burst_detected = true
        end
    end

    if !burst_detected
        println("✅ Degradation remains within safe biological limits.")
    end

    # Create results directory if it does not exist
    results_dir = joinpath(@__DIR__, "..", "results")
    mkpath(results_dir)

    # Generate the clinical threshold plot
    p = plot(time_days, zn_concentration_proxy, 
             title="Zn²⁺ Interface Accumulation vs Toxicity Limit",
             xlabel="Time (days)", 
             ylabel="Local Ion Concentration Proxy",
             color=:red, lw=2, label="Ion Accumulation")

    hline!(p, [toxicity_limit], ls=:dash, color=:black, lw=2, label="Biological Safety Limit")

    output_path = joinpath(results_dir, "toxicity_burst_analysis.png")
    savefig(p, output_path)

    println("\n✅ Toxicity analysis complete! Artifact saved as '", output_path, "'")
end

# Execute the simulation
run_toxicity_analysis()