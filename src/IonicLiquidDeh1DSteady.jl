module IonicLiquidDeh1DSteady
    using PrecompileTools: @setup_workload, @compile_workload  
    using StaticArrays
    using CoolProp
    using Interpolations
    using NonlinearSolve
    using BoundaryValueDiffEq
    using UnPack: @unpack

    include("TypesInterface/types_interface.jl")
    include("FluidProperties/Props.jl")
    include("SolverX/solver.jl")
    include("Simulation/simulation.jl")

    @setup_workload begin
        @compile_workload begin
            include("Simulation/_precompilation.jl")
        end
    end


end


