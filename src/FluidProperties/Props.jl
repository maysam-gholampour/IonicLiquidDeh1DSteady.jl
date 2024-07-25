export IonicLiquid, CreCOPlus5100

abstract type IonicLiquid end
struct CreCOPlus5100 <: IonicLiquid end


include("Air.jl")
include("Water.jl")
include("CreCOPlus5100.jl")


