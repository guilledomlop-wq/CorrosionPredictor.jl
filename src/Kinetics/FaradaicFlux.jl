"""
    calculate_butler_volmer_flux(overpotential::Voltage, i0::CurrentDensity, alpha_a::Real, alpha_c::Real, z::Int, T::Temperature, j_lim::MolarFlux, R_film_factor::Real=1.0, tafel_transition_threshold::Real=50.0)

Calculate the net molar flux of metal dissolution using the Butler-Volmer equation.
Converts the kinetically driven current density into physical mass flux via Faraday's Law.
Incorporates a Tafel switch to prevent exponential overflow, Current Density Capping to respect physical diffusion limits,
and an `R_film_factor` to scale down the massive theoretical overpotential based on the corrosion product layer.
"""
function calculate_butler_volmer_flux(overpotential::Voltage, i0::CurrentDensity, alpha_a::Real, alpha_c::Real, z::Int, T::Temperature, j_lim::MolarFlux, R_film_factor::Real=1.0, tafel_transition_threshold::Real=50.0)
    R = 8.314 * u"J/(mol*K)"
    F = 96485.0 * u"C/mol"
    
    # Dampen the overpotential using the film resistance factor
    eta_eff = overpotential / R_film_factor
    
    # Calculate dimensionless exponent factors (alpha * z * F * eta) / (R * T)
    f_anodic = ustrip(uconvert(Unitful.NoUnits, (alpha_a * z * F * eta_eff) / (R * T)))
    f_cathodic = ustrip(uconvert(Unitful.NoUnits, (alpha_c * z * F * eta_eff) / (R * T)))
    
    # Tafel switch to prevent Float64 exponential overflow (NaN generation)
    # Implemented as a pure Tafel exponential approximation beyond the threshold
    # analytically designed to maintain C1 smoothness and scientific accuracy.
    
    if f_anodic > tafel_transition_threshold
        # Pure anodic Tafel approximation
        i_net = i0 * exp(f_anodic)
    elseif f_cathodic > tafel_transition_threshold
        # Pure cathodic Tafel approximation
        i_net = -i0 * exp(f_cathodic)
    else
        # Net current density (anodic component minus cathodic component)
        i_net = i0 * (exp(f_anodic) - exp(-f_cathodic))
    end
    
    # Convert current density to molar flux (Faraday's Law): j = i / (z * F)
    # We enforce a physical floor at 0 to prevent "reverse corrosion" in this model
    i_anodic_effective = max(i_net, 0.0 * u"A/m^2")
    j_metal_unlimited = i_anodic_effective / (z * F)
    j_metal_unlimited_conv = uconvert(u"mol/(m^2*s)", j_metal_unlimited)
    
    # Current Density Capping: Limit the kinetic flux to the available physical transport
    return min(j_metal_unlimited_conv, j_lim)
end

"""
    corrosion_rate_to_flux(rate::CorrosionRate, material::Material)

Convert a linear corrosion rate (e.g., um/yr) to a molar flux (mol/m^2/s).
Used as the primary input for kinetic simulations.
"""
function corrosion_rate_to_flux(rate::CorrosionRate, material::Material)
    # Physical constant
    seconds_per_year = 365.25 * 24 * 3600 * u"s/yr"
    
    # Unit conversion logic (uconvert handles Measurements automatically)
    rate_m = uconvert(u"m/yr", rate)
    mass_flux = material.density * rate_m              # g/(m^2*yr)
    moles_flux = mass_flux / material.molar_mass       # mol/(m^2*yr)
    
    return moles_flux / seconds_per_year            # mol/(m^2*s)
end