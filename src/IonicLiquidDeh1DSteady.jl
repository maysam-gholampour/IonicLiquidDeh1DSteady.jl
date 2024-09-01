module IonicLiquidDeh1DSteady
    using StaticArrays
    using CoolProp
    using Interpolations
    using NonlinearSolve
    using BoundaryValueDiffEq
    using SciMLBase

    include("TypesInterface/types_interface.jl")
    include("FluidProperties/Props.jl")
    include("SolverX/solver.jl")
    include("Simulation/simulation.jl")

end
