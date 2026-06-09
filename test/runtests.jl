# =========================================================================
# test/runtests.jl
# Comprehensive test suite for the CorrosionPredictor core mechanics.
# =========================================================================

using Test
using CorrosionPredictor
using Unitful

@testset "CorrosionPredictor.jl Test Suite" begin

    @testset "1. Materials and Type System" begin
        # Verify strict typing and immutable constants
        @test ZINC.valency == 2
        @test ZINC.name == "Zinc"
        @test MAGNESIUM.valency == 2
        
        # Verify dimension consistency
        @test typeof(ZINC.standard_potential) <: Voltage
        @test typeof(ZINC.density) <: Density
    end

    @testset "2. Physics Core" begin
        # Test dimensionless film resistance scaling
        # At 1 micrometer (1e-6 m), the resistance factor should be base (15.0) + 1.0 = 16.0
        R_film = calculate_film_resistance_factor(1e-6u"m")
        @test isapprox(R_film, 16.0, atol=1e-5)
        
        # Test interface concentration physical boundaries (cannot drop below zero)
        c_bulk = 10.0u"mol/m^3"
        D = 1e-9u"m^2/s"
        delta = 10e-6u"m"
        
        # Massive consumption flux that would theoretically cause negative concentration
        extreme_flux = -100.0u"mol/(m^2*s)" 
        c_int = calculate_interface_concentration(c_bulk, extreme_flux, D, delta)
        
        @test c_int == 0.0u"mol/m^3"
    end

    @testset "3. Kinetics Submodule" begin
        # Test basic thermodynamics: Debye-Hückel / Davies activity coefficient
        gamma_z2 = calculate_activity_coefficient(0.15, 2)
        @test gamma_z2 > 0.0 && gamma_z2 < 1.0 # Should reduce activity in physiological saline
        
        # Test Faraday conversion logic (Mass flux vs Current density)
        rate = 20.0u"μm/yr"
        flux = corrosion_rate_to_flux(rate, ZINC)
        @test typeof(flux) <: MolarFlux
        @test ustrip(flux) > 0.0
    end

end