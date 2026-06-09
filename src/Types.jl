# =========================================================================
# src/Types.jl
# Centralized Type System and Unitful Dimensions for the Predictive Corrosion Platform.
# Eliminates duplication and ensures strict dimensional consistency for SciML solvers.
# =========================================================================

using Unitful

# =========================================================================
# CENTRALIZED PHYSICAL DIMENSIONS
# Using dimension() to ensure Type Safety without Unicode warnings in IDEs
# =========================================================================

# Base Physical and Material Dimensions
const LENGTH_D = Unitful.dimension(u"m")
const AREA_D = Unitful.dimension(u"m^2")
const VOLUME_D = Unitful.dimension(u"m^3")
const DENSITY_D = Unitful.dimension(u"g/cm^3")
const MOLAR_MASS_D = Unitful.dimension(u"g/mol")

# Transport and Concentration Dimensions
const CONCENTRATION_D = Unitful.dimension(u"mol/m^3")
const DIFFUSION_D = Unitful.dimension(u"m^2/s")
const MOLAR_FLUX_D = Unitful.dimension(u"mol/(m^2*s)")

# Kinetics and Degradation Dimensions
const CORROSION_RATE_D = Unitful.dimension(u"m/s")
const MASS_RATE_D = Unitful.dimension(u"g/s")

# Electrochemical Dimensions
const VOLTAGE_D = Unitful.dimension(u"V")
const CURRENT_DENSITY_D = Unitful.dimension(u"A/m^2")
const GAS_CONSTANT_D = Unitful.dimension(u"J/(mol*K)")
const FARADAY_CONSTANT_D = Unitful.dimension(u"C/mol")
const TEMPERATURE_D = Unitful.dimension(u"K")
const CAPACITANCE_PER_AREA_D = Unitful.dimension(u"F/m^2")

# =========================================================================
# STRICT TYPE ALIASES
# Type aliases for strict dimension validation across SciML and PDE functions
# =========================================================================

const Length{T} = Quantity{T, LENGTH_D}
const Area{T} = Quantity{T, AREA_D}
const FluidVolume{T} = Quantity{T, VOLUME_D}
const Density{T} = Quantity{T, DENSITY_D}
const MolarMass{T} = Quantity{T, MOLAR_MASS_D}

const Concentration{T} = Quantity{T, CONCENTRATION_D}
const DiffusionCoefficient{T} = Quantity{T, DIFFUSION_D}
const MolarFlux{T} = Quantity{T, MOLAR_FLUX_D}

const CorrosionRate{T} = Quantity{T, CORROSION_RATE_D}
const MassRate{T} = Quantity{T, MASS_RATE_D}

const Voltage{T} = Quantity{T, VOLTAGE_D}
const CurrentDensity{T} = Quantity{T, CURRENT_DENSITY_D}
const GasConstant{T} = Quantity{T, GAS_CONSTANT_D}
const FaradayConstant{T} = Quantity{T, FARADAY_CONSTANT_D}
const Temperature{T} = Quantity{T, TEMPERATURE_D}
const CapacitancePerArea{T} = Quantity{T, CAPACITANCE_PER_AREA_D}