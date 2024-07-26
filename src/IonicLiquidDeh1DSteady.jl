module IonicLiquidDeh1DSteady
    using StaticArrays
    using CoolProp
    using Interpolations
    using NonlinearSolve
    using BoundaryValueDiffEq

    include("FluidProperties/Props.jl")
    include("solver.jl")

end
