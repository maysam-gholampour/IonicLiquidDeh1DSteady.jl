#= 
    Manuel R. Conde,
    Properties of aqueous solutions of lithium and calcium chlorides: formulations for use in air conditioning equipment design,
    International Journal of Thermal Sciences,
    Volume 43, Issue 4,
    2004,
    Pages 367-382,
    ISSN 1290-0729,
    https://doi.org/10.1016/j.ijthermalsci.2003.09.003.
    (https://www.sciencedirect.com/science/article/pii/S1290072903001625)
    Abstract: The dehydration of air, for air conditioning purposes, either for human comfort or for industrial processes, is done most of the times by making it contact a surface at a temperature below its dew point. In this process not only is it necessary to cool that surface continuously, but also the air is cooled beyond the temperature necessary to the process, thus requiring reheating after dehumidification. Although the equipment for this purpose is standard and mostly low-cost, the running costs are high and high grade energy is dissipated at very low efficiency. Alternative sorption-based processes require only low grade energy for regeneration of the sorbent materials, thus incurring lower running costs. On the other hand, sorption technology equipment is usually more expensive than standard mechanical refrigeration equipment, which is essentially due to their too small market share. This paper reports the development of calculation models for the thermophysical properties of aqueous solutions of the chlorides of lithium and calcium, particularly suited for use as desiccants in sorption-based air conditioning equipment. This development has been undertaken in order to create consistent methods suitable for use in the industrial design of liquid desiccant-based air conditioning equipment. We have reviewed sources of measured data from 1850 onwards, and propose calculation models for the following properties of those aqueous solutions: Solubility boundary, vapour pressure, density, surface tension, dynamic viscosity, thermal conductivity, specific thermal capacity and differential enthalpy of dilution.
    Keywords: Liquid desiccants; Properties; Air conditioning; Open absorption; Lithium chloride; Calcium chloride; Calculation models
=#

@testset "LiCl Properties" begin
    T = 50.0 + 273.15
    ξ = 0.3
    # Fig. 3. Relative vapour pressure of aqueous solutions of lithium chloride
    @test _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ, LiCl()) / P_H2O(T) ≈ 0.44593378 atol=1e-8

    # Fig. 6. Relative densities of aqueous solutions of lithium chloride.
    @test _ρₛₒₗ(T, ξ, LiCl()) / ρ_H2O(1 - T / 647.226)  ≈ 1.18397 atol=1e-5

    # Fig. 10. Dynamic viscosity of aqueous solutions of lithium chloride.
    @test _μₛₒₗ(T, ξ, LiCl()) ≈ 0.002049 atol=1e-6

    # Fig. 17 and 18. Specific thermal capacity of aqueous solutions of lithium chloride
    @test _cpₛₒₗ(T, ξ, LiCl()) ≈ 3016.2847 atol=1e-4

    # Fig. 13. Thermal conductivity of aqueous solutions of calcium chloride.
    @test _𝑘ₛₒₗ(T, ξ, LiCl()) ≈ 0.580522 atol=1e-6

    # Fig. 8. Relative values of the surface tension of aqueous solutions of lithium chloride
    @test _σₛₒₗ(T, ξ, LiCl()) / σ_H2O(T / 647.226) ≈ 1.266139 atol=1e-6

    # Fig. 21. Calculated and measured values of the differential enthalpy of
    # dilution for aqueous solutions of lithium chloride.
    @test _iₛₒₗ(T, ξ, LiCl()) ≈ 116705.086 atol=1e-3


    @test calculate_T_sol(116705.086, ξ, LiCl()) ≈ 50.0 + 273.15 atol=1e-3
end
