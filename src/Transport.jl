# =========================================================================
# src/Transport.jl
# =========================================================================

using Unitful
using ModelingToolkit, MethodOfLines
using DomainSets

"""
    calculate_diffusive_flux(c_bulk::Concentration, c_interface::Concentration, D::DiffusionCoefficient, delta::Length)

Calculate the diffusive flux based on the concentration gradient: j = D * (c_bulk - c_interface) / delta.
"""
function calculate_diffusive_flux(c_bulk::Concentration, c_interface::Concentration, D::DiffusionCoefficient, delta::Length)
    return D * (c_bulk - c_interface) / delta
end

"""
    build_diffusion_model(D::DiffusionCoefficient, delta::Length, c_bulk::Concentration)

Initializes the symbolic components for a 1D spatio-temporal diffusion PDE system.
Prepares the variables (t, x, c) and the core equation for MethodOfLines.
"""
function build_diffusion_model(D::DiffusionCoefficient, delta::Length, c_bulk::Concentration)
    @parameters t x
    @variables c(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2

    # PDE: ∂c/∂t = D * ∂²c/∂x²
    eq = Dt(c(t, x)) ~ ustrip(uconvert(u"m^2/s", D)) * Dxx(c(t, x))

    # Spatial domain: x ∈ [0, delta]
    domain = [t ∈ Interval(0.0, Inf),
              x ∈ Interval(0.0, ustrip(uconvert(u"m", delta)))]

    return eq, domain, t, x, c
end