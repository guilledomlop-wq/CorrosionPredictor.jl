# =========================================================================
# src/Kinetics/MixedPotentialSolver.jl
# Component sub-file for Kinetics root-finding and mixed potential tracking.
# Modified to ensure electrochemical keys override static corrosion fallbacks.
# =========================================================================

"""
    get_dynamic_o2_potential(pH_local::Real, T_K::Real, E_eq_base_74::Real)

Calculates the shift in the Oxygen reduction equilibrium potential 
due to the accumulation of OH- ions (alkalization) at the stent interface.
"""
function get_dynamic_o2_potential(pH_local::Real, T_K::Real, E_eq_base_74::Real)
    R = 8.314
    F = 96485.0
    
    # Nernst slope for OH- generation: ~61.5 mV drop per pH unit at 37 C
    nernst_slope = (2.303 * R * T_K) / F
    
    # Calculate the shift from the physiological 7.4 baseline
    pH_shift = pH_local - 7.4
    
    return E_eq_base_74 - (nernst_slope * pH_shift)
end

"""
    get_dynamic_diffusion_limit(pH_local::Real, i_lim_base::Real)

Simulates the physical blockage of the stent surface due to the 
precipitation of Zn(OH)2 and ZnO in highly alkaline environments.
Reduces the available diffusion limit using a sigmoidal decay.
"""
function get_dynamic_diffusion_limit(pH_local::Real, i_lim_base::Real)
    # Critical pH where massive precipitation forms a dense passivation layer
    pH_crit = 9.0 
    
    # Steepness of the precipitation transition
    k_precip = 3.0 
    
    # Sigmoidal drop simulating area blockage: Drops to 5% of original capacity at high pH
    # Floor of 0.05 prevents boundary conditions from dropping to or below zero
    block_factor = 0.05 + 0.95 / (1.0 + exp(k_precip * (pH_local - pH_crit)))
    
    return i_lim_base * block_factor
end

"""
    find_mixed_potential(E_eq_a::Voltage, i0_a::CurrentDensity, alpha_a::Real, z_a::Int, j_lim_a::MolarFlux, E_eq_c_base::Voltage, c_bulk_c::Concentration, i0_c::CurrentDensity, alpha_c::Real, f_2e::Real, T::Temperature, R_film_factor::Real=1.0, pH_local::Real=7.4)

Finds the dynamic Mixed Potential (E_corr) using Roots.jl, where the sum of anodic and cathodic currents is zero.
Enforces a strict physical ceiling based on the diffusion limit.
Now calculates a dynamic, concentration-coupled Nernst potential for the cathodic reaction (oxygen reduction) internally.
Introduces `R_film_factor` to scale the effective overpotential based on the protective oxide layer thickness.
Incorporates dynamic pH tracking to shift the O2 potential and physical precipitation blocking to choke the diffusion limit.
[NEW] Dynamically couples the cathodic valency to the active material valency and the f_2e fraction.
"""
function find_mixed_potential(E_eq_a::Voltage, i0_a::CurrentDensity, alpha_a::Real, z_a::Int, j_lim_a::MolarFlux,
                              E_eq_c_base::Voltage, c_bulk_c::Concentration, i0_c::CurrentDensity, alpha_c::Real, f_2e::Real, T::Temperature, R_film_factor::Real=1.0, pH_local::Real=7.4)
    R = 8.314
    F = 96485.0

    # Dynamic effective cathodic valency based on the 2-electron pathway fraction
    z_c_effective = 4.0 * (1.0 - f_2e) + 2.0 * f_2e

    # Strip units for internal root-finding optimizations
    T_num = ustrip(uconvert(u"K", T)) # Standardizing variables via uconvert to absolute temperature scale
    E_eq_a_num = ustrip(uconvert(u"V", E_eq_a))
    i0_a_num = ustrip(uconvert(u"A/m^2", i0_a))
    j_lim_a_num = ustrip(uconvert(u"mol/(m^2*s)", j_lim_a))
    i_lim_a_num = j_lim_a_num * (z_a * F) # Faraday conversion to current density ceiling
    
    # Apply dynamic pH physical precipitation blockage to the anodic diffusion limit
    i_lim_a_dynamic = get_dynamic_diffusion_limit(pH_local, i_lim_a_num)
    
    E_eq_c_base_num = ustrip(uconvert(u"V", E_eq_c_base))
    c_bulk_c_num = ustrip(uconvert(u"mol/m^3", c_bulk_c))
    i0_c_num = ustrip(uconvert(u"A/m^2", i0_c))
    
    # Standard state reference for thermodynamic normalization
    c0_num = 1000.0 # Standard concentration 1M = 1000 mol/m^3
    ionic_strength_M = 0.15 # Physiological ionic strength
    
    # Calculate activity coefficient for the neutral cathodic species (O2, z=0)
    gamma_c = calculate_activity_coefficient(ionic_strength_M, 0)
    
    # Anodic current (with Tafel switch, strict capping, and film resistance damping)
    function anodic_current(E)
        eta = E - E_eq_a_num
        eta_eff = eta / R_film_factor # Apply physical damping from the corrosion product layer
        f_anodic = (alpha_a * z_a * F * eta_eff) / (R * T_num)
        threshold = 50.0 
        
        if f_anodic > threshold
            # First-order Taylor expansion at f = threshold to prevent evaluation crash
            i_kinetic = i0_a_num * exp(threshold) * (1.0 + (f_anodic - threshold))
        else
            i_kinetic = i0_a_num * (exp(f_anodic) - exp(-f_anodic))
        end
        
        i_effective = max(i_kinetic, 0.0)
        return min(i_effective, i_lim_a_dynamic)
    end
    
    # Cathodic current (Oxygen reduction with dynamic Nernst potential and film resistance damping)
    function cathodic_current(E, E_eq_c_dynamic)
        eta = E - E_eq_c_dynamic
        eta_eff = eta / R_film_factor # Apply physical damping from the corrosion product layer
        
        # Define full Butler-Volmer exponent factors to mathematically mirror the anodic shielding
        f_cathodic = (alpha_c * z_c_effective * F * eta_eff) / (R * T_num)
        f_anodic = ((1.0 - alpha_c) * z_c_effective * F * eta_eff) / (R * T_num)
        
        threshold = 50.0
        
        if f_anodic > threshold
            # Extreme anodic polarization: C1 continuous linear extrapolation
            i_kinetic = i0_c_num * exp(threshold) * (1.0 + (f_anodic - threshold))
        elseif -f_cathodic > threshold
            # Extreme cathodic polarization: C1 continuous linear extrapolation
            i_kinetic = -i0_c_num * exp(threshold) * (1.0 + (-f_cathodic - threshold))
        else
            # Safe equilibrium zone: Full Butler-Volmer evaluates perfectly
            i_kinetic = i0_c_num * (exp(f_anodic) - exp(-f_cathodic))
        end
        
        return i_kinetic
    end
    
    # Objective function: Find voltage E where total current is zero
    function total_current(E)
        i_a = anodic_current(E)
        
        # Dynamic coupling: Convert anodic current back to metal flux to determine oxygen depletion
        j_metal_num = i_a / (z_a * F)
        
        # Reconstruct stoichiometric fraction dynamically matching the active material parameters
        # j_lim_a_num maps the metal flux equivalent supported by the available O2 diffusion boundary layer limit
        j_O2_limit_num = j_lim_a_num * ((z_a / 4) * (1.0 - f_2e) + (z_a / 2) * f_2e) # Equivalent O2 baseline limit reference
        j_O2_consumption_num = j_metal_num * ((z_a / 4) * (1.0 - f_2e) + (z_a / 2) * f_2e) # Direct stoichiometric consumption flux
        
        # Calculate dynamic interface concentration of the cathodic species (O2)
        # Using the linear diffusion limit relationship: c_interface = c_bulk * (1 - j_O2 / j_O2_limit)
        # We enforce a biological trace floor to prevent log(0) singularity
        c_interface_c_num = max(c_bulk_c_num * (1.0 - (j_O2_consumption_num / j_O2_limit_num)), 1e-6)
        
        # Dynamically calculate the baseline O2 Nernst potential based on local pH
        E_eq_c_base_dynamic = get_dynamic_o2_potential(pH_local, T_num, E_eq_c_base_num)
        
        # Dynamic Nernst potential for the cathodic reaction based on local depletion
        # Scaled dynamically with Davies activity and standard state references
        activity_ratio_c = (c_interface_c_num / c0_num) * gamma_c
        E_eq_c_dynamic = E_eq_c_base_dynamic + (R * T_num / (z_c_effective * F)) * log(activity_ratio_c)
        
        i_c = cathodic_current(E, E_eq_c_dynamic)
        return i_a + i_c
    end
    
    # Determine maximum possible cathodic equilibrium potential (at zero current, when interface = bulk)
    # This serves as a mathematically absolute and safe upper bound for the Bisection method
    activity_ratio_c_max = (max(c_bulk_c_num, 1e-6) / c0_num) * gamma_c
    E_eq_c_base_dynamic_max = get_dynamic_o2_potential(pH_local, T_num, E_eq_c_base_num)
    E_eq_c_max = E_eq_c_base_dynamic_max + (R * T_num / (z_c_effective * F)) * log(activity_ratio_c_max)
    
    # Execute Bisection root-finding within the thermodynamic limits
    E_corr_num = find_zero(total_current, (E_eq_a_num, E_eq_c_max), Bisection())
    
    return E_corr_num * u"V"
end

"""
    calculate_stoichiometric_fluxes(j_metal::MolarFlux, f_2e::Real, material::Material{T}) where T <: BioresorbableMaterial

Calculate the stoichiometric fluxes for oxygen consumption and hydrogen peroxide generation
based on the electrochemical kinetics of metal corrosion. 
"""
function calculate_stoichiometric_fluxes(j_metal::MolarFlux, f_2e::Real, material::Material{T}) where T <: BioresorbableMaterial
    valency = material.valency

    if !(0.0 <= f_2e <= 1.0)
        throw(ArgumentError("f_2e must be between 0.0 and 1.0"))
    end

    j_metal_conv = uconvert(u"mol/m^2/s", j_metal)

    O2_consumption = ((valency / 4) * (1 - f_2e) + (valency / 2) * f_2e) * j_metal_conv
    H2O2_generation = (valency / 2) * f_2e * j_metal_conv

    return (O2_consumption = O2_consumption, H2O2_generation = H2O2_generation)
end

"""
    calculate_surface_environment(material::Material, corrosion_rate::CorrosionRate, f_2e::Real, env_params::Dict)

Integrate kinetics and transport logic to determine the environment at the stent interface.
Calculates interface concentrations (c_interface) for O2, Zn2+, and H2O2.
Supports both static reference corrosion rates and dynamic Mixed Potential kinetics.
"""
function calculate_surface_environment(material::Material, corrosion_rate::CorrosionRate, f_2e::Real, env_params::Dict)
    delta = env_params["boundary_layer_thickness"]
    c_O2_bulk = env_params["c_O2_bulk"]
    D_O2 = env_params["D_O2"]
    
    j_O2_limit = (D_O2 * c_O2_bulk) / delta
    
    valency = material.valency
    stoich_factor = (valency / 4) * (1 - f_2e) + (valency / 2) * f_2e
    j_metal_limit = j_O2_limit / stoich_factor

    if haskey(env_params, "i0_zinc") || haskey(env_params, "i0_magnesium")
        T = env_params["T_temperature"]
        i0_metal = material.name == "Zinc" ? env_params["i0_zinc"] : env_params["i0_magnesium"]
        alpha_a = env_params["alpha_a"]
        alpha_c = env_params["alpha_c"]
        
        current_L = get(env_params, "current_layer_thickness", 0.0u"m")
        R_film_factor = calculate_film_resistance_factor(current_L)
        pH_local = get(env_params, "pH_local", 7.4)
        
        E_eq_o2_base = 0.40u"V"
        i0_o2 = 1e-4u"A/m^2"
        
        ionic_strength = get(env_params, "ionic_strength_M", 0.15)
        c_metal_bulk = env_params["c_Zn2_bulk"]
        E_eq_metal = calculate_nernst_potential(material, c_metal_bulk, T, ionic_strength)
        
        E_corr = find_mixed_potential(E_eq_metal, i0_metal, alpha_a, valency, j_metal_limit,
                                      E_eq_o2_base, c_O2_bulk, i0_o2, alpha_c, f_2e, T, R_film_factor, pH_local)
                                      
        eta = E_corr - E_eq_metal
        j_metal = calculate_butler_volmer_flux(eta, i0_metal, alpha_a, alpha_c, valency, T, j_metal_limit, R_film_factor)
    else
        j_metal_static = corrosion_rate_to_flux(corrosion_rate, material)
        j_metal = min(j_metal_static, j_metal_limit)
    end
    
    fluxes = calculate_stoichiometric_fluxes(j_metal, f_2e, material)
    
    c_metal_bulk = env_params["c_Zn2_bulk"] 
    D_metal = env_params["D_Zn2"]
    c_metal_interface = calculate_interface_concentration(c_metal_bulk, j_metal, D_metal, delta)
    
    c_O2_interface = calculate_interface_concentration(c_O2_bulk, -fluxes.O2_consumption, D_O2, delta)
    
    c_H2O2_bulk = env_params["c_H2O2_bulk"]
    D_H2O2 = env_params["D_H2O2"]
    c_H2O2_interface = calculate_interface_concentration(c_H2O2_bulk, fluxes.H2O2_generation, D_H2O2, delta)
    
    return (
        j_metal = j_metal,
        fluxes = fluxes,
        c_Zn2_interface = c_metal_interface,
        c_O2_interface = c_O2_interface,
        c_H2O2_interface = c_H2O2_interface
    )
end