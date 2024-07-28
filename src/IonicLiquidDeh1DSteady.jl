module IonicLiquidDeh1DSteady
    using StaticArrays
    using CoolProp
    using Interpolations
    using NonlinearSolve
    using BoundaryValueDiffEq

    include("FluidProperties/Props.jl")
    include("SolverX/solver.jl")

end
