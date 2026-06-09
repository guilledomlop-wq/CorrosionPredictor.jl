using Unitful
using Roots # Explicit dependency for root finding algorithms

"""
    calculate_activity_coefficient(ionic_strength_M::Real, z::Int)

Calculates the activity coefficient (gamma) using the Davies equation,
which extends Debye-Hückel for physiological ionic strengths (~0.15 M).
"""
function calculate_activity_coefficient(ionic_strength_M::Real, z::Int)
    A = 0.509 # Constant for water at 37 C
    I = ionic_strength_M
    
    # Davies equation
    log_gamma = -A * z^2 * (sqrt(I) / (1.0 + sqrt(I)) - 0.3 * I)
    return 10.0^log_gamma
end

"""
    calculate_nernst_potential(material::Material, c_ion_interface::Concentration, T::Temperature, ionic_strength_M::Real=0.15)

Calculate the dynamic equilibrium potential (Nernst potential) for the metal dissolution reaction.
Adjusts the standard reduction potential based on local metal ion accumulation and real chemical activity.
"""
function calculate_nernst_potential(material::Material, c_ion_interface::Concentration, T::Temperature, ionic_strength_M::Real=0.15)
    # Universal constants hardcoded for internal consistency if not passed
    R = 8.314 * u"J/(mol*K)"
    F = 96485.0 * u"C/mol"
    z = material.valency
    E0 = material.standard_potential
    
    # Standard reference concentration (1 Molar) to ensure dimensionless log argument
    c0 = 1000.0 * u"mol/m^3"
    
    # Avoid log(0) domain error by applying a trace background concentration if exactly 0
    safe_c = c_ion_interface <= 0.0u"mol/m^3" ? 1e-9u"mol/m^3" : c_ion_interface
    
    # Calculate activity coefficient based on Debye-Hückel / Davies
    gamma = calculate_activity_coefficient(ionic_strength_M, z)
    
    # Activity ratio (concentration ratio * activity coefficient for physiological solutions)
    activity_ratio = ustrip(uconvert(Unitful.NoUnits, safe_c / c0)) * gamma
    
    # Nernst Equation: E_eq = E0 + (RT/zF) * ln(a_M)
    E_eq = E0 + (R * T / (z * F)) * log(activity_ratio)
    
    return uconvert(u"V", E_eq)
end