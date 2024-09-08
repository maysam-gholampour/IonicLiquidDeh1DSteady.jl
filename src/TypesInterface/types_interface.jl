export IonicLiquid, CreCOPlus5100, LiCl, CaCl2
export Coil, PlateFinCircularTube
export InputData, FluidThermalData

# ==================== IonicLiquid ====================
abstract type IonicLiquid end
struct CreCOPlus5100 <: IonicLiquid end
struct LiCl <: IonicLiquid end
struct CaCl2<:IonicLiquid end
# ==================== Coil ====================
abstract type Coil end
struct PlateFinCircularTube{T1<:AbstractFloat, T2<:Integer} <: Coil 
    δ_fin::T1
    D_tube_outside::T1
    N_tube_per_row::T2
    N_row::T2
    N_fin::T2
    FD::T1
    H::T1
    FS::T1
    σ::T1
end
# ==================== InputData ====================
abstract type InputData end
struct FluidThermalData{T1<:AbstractFloat,T2<:IonicLiquid}<:InputData
    T_air::T1
    T_wb_air::T1
    m_dot_air::T1
    m_dot_sol::T1
    T_sol_in::T1
    X_sol_in::T1
    T_ref_in::T1  #DELME: delete them
    T_ref_out::T1  #DELME: delete them
    IL::T2
    Q::T1
    P_sat::T1  #DELME: delete them
    T_sat::T1  #DELME: delete them
    Le::T1
end