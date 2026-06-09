# =========================================================================
# src/PhysicsCore.jl
# Low-level physical and geometric helper functions.
# Evaluated before Kinetics logic to prevent circular dependencies.
# =========================================================================

using Unitful

"""
    calculate_film_resistance_factor(L)

Calculates the dimensionless film resistance damping factor based on the current
corrosion product layer thickness (L). This factor represents the ohmic and physical
damping from the protective oxide layer.
"""
function calculate_film_resistance_factor(L)
    return 15.0 + ustrip(uconvert(Unitful.NoUnits, L / 1e-6u"m"))
end

"""
    calculate_interface_concentration(c_bulk, j_surface, D, delta)

Calculate the steady-state concentration at the material-fluid interface (c_interface) 
based on Fick's first law: c_interface = c_bulk + (j_surface * delta) / D.
"""
function calculate_interface_concentration(c_bulk, j_surface, D, delta)
    # Core physics: c_int = c_bulk + (flux * distance) / Diffusion 
    c_int = c_bulk + (j_surface * delta) / D
    
    # Physical constraint: concentration cannot be negative 
    # Type safety ensures c_int and zero_val share compatible units automatically
    zero_val = 0.0 * unit(c_int)
    return c_int < zero_val ? zero_val : c_int
end