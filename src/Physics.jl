# =========================================================================
# src/Physics.jl
# =========================================================================

using Unitful
using Unitful: g, cm, mol, s, μm, hr, mg, d
using Unitful.DefaultSymbols

"""
    calculate_mass_loss_rate(j_metal::MolarFlux, area::Area, molar_mass::MolarMass)

Calculate the mass change rate (dM/dt) of a metal due to corrosion.

# Arguments
- `j_metal::MolarFlux`: Molar flux of the metal ions (mol/(m^2*s)).
- `area::Area`: Surface area of the corroding metal (m^2).
- `molar_mass::MolarMass`: Molar mass of the metal (g/mol).

# Returns
- `MassRate`: The rate of mass change (g/s).
"""
function calculate_mass_loss_rate(j_metal::MolarFlux, area::Area, molar_mass::MolarMass)
    return (j_metal * area * molar_mass) |> u"g/s"
end

"""
    calculate_layer_growth_rate(dm_dt::MassRate, k_growth::Real, density::Density, area::Area)

Calculate the thickness change rate (dL/dt) of a corrosion product layer.

# Arguments
- `dm_dt::MassRate`: The rate of mass change (e.g., g/s), typically mass loss of the metal.
- `k_growth::Real`: A dimensionless constant representing the volumetric ratio of the corrosion product to the consumed metal.
- `density::Density`: Density of the corrosion product (g/cm^3).
- `area::Area`: Surface area where the layer is growing (m^2).

# Returns
- `CorrosionRate`: The rate of thickness change (e.g., μm/year).
"""
function calculate_layer_growth_rate(dm_dt::MassRate, k_growth::Real, density::Density, area::Area)
    return (k_growth * dm_dt / density / area) |> u"μm/yr"
end