using Unitful

"""
    BioresorbableMaterial

Abstract type representing the category of biodegradable metals for medical applications.
Used to facilitate multiple dispatch in kinetic and diffusion simulations.
"""
abstract type BioresorbableMaterial end

"""
    ZincType <: BioresorbableMaterial
    MagnesiumType <: BioresorbableMaterial

Concrete marker types for specific material identities.
"""
struct ZincType <: BioresorbableMaterial end
struct MagnesiumType <: BioresorbableMaterial end

"""
    Material{T <: BioresorbableMaterial}

Represents a material with physical properties, carrying its type identity `T` 
for dispatch purposes.

# Fields
- `name::String`: The name of the material
- `density::Density`: Density (strictly validated dimensions)
- `molar_mass::MolarMass`: Molar mass (strictly validated dimensions)
- `valency::Int`: Valency (number of electrons involved in oxidation)
- `standard_potential::Voltage`: Standard reduction potential E^0 vs SHE (strictly validated dimensions)
- `alkalization_slope::Float64`: Stoichiometric accumulation factor for local pH
- `max_ph_limit::Float64`: Absolute saturation limit for local alkalization
- `base_oxygen_potential::Voltage`: Baseline O2 reduction potential at standard physiological pH
- `film_resistance_base::Float64`: Base dimensionless ohmic damping from the pristine oxide layer
"""
struct Material{T <: BioresorbableMaterial}
    name::String
    density::Density
    molar_mass::MolarMass
    valency::Int
    standard_potential::Voltage
    alkalization_slope::Float64
    max_ph_limit::Float64
    base_oxygen_potential::Voltage
    film_resistance_base::Float64
end

# Material constants
const ZINC = Material{ZincType}(
    "Zinc",
    7.14u"g/cm^3",
    65.38u"g/mol",
    2,
    -0.7618u"V", # Standard potential for Zn2+ + 2e- -> Zn
    0.15,        # alkalization_slope
    10.5,        # max_ph_limit
    0.40u"V",    # base_oxygen_potential
    15.0         # film_resistance_base
)

const MAGNESIUM = Material{MagnesiumType}(
    "Magnesium",
    1.74u"g/cm^3",
    24.305u"g/mol",
    2,
    -2.372u"V",  # Standard potential for Mg2+ + 2e- -> Mg
    0.30,        # alkalization_slope (Estimated stoichiometric difference)
    11.5,        # max_ph_limit (Higher alkaline saturation for Mg)
    0.40u"V",    # base_oxygen_potential
    10.0         # film_resistance_base (Different native oxide properties)
)

export Material, BioresorbableMaterial, ZincType, MagnesiumType, ZINC, MAGNESIUM