module IonicLiquidDeh1DSteady
    __precompile__(true)
    using Reexport: @reexport
    
    using PrecompileTools: @setup_workload, @compile_workload  
    using StaticArrays
    @reexport using CoolProp
    @reexport using NonlinearSolve
    @reexport using UnPack: @unpack
    using Interpolations
    using BoundaryValueDiffEq

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


