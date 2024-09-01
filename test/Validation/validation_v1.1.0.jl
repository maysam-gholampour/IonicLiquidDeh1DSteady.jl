using IonicLiquidDeh1DSteady
using CoolProp


begin 
    δ_fin_evap = 0.00013
    D_tube_outside_evap = 0.0101
    N_tube_per_row_evap = 6
    N_row_evap = 6
    N_fin_evap = 48
    FD_evap = 0.155
    H_evap = 0.132
    FS_evap = 0.00254
    Le = 0.86
    σ_evap = 1.0
    # ========================================
    IL = CreCOPlus5100()
    T_air_amb = 25.48 + 273.15 + 1.5 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 21.21 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 1100 * 2.57 / 60 / 1000  # kg/s  #FIXME: 1100 is density
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.76
    T_ref_in_evap = 12.19 + 273.15  # K
    T_ref_out_evap = 17.04 + 273.15  # K
    Q_evap = 1192.0
    P_evap = 4.3e5
    T_sat_evap = 10.51 + 273.15
end

plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
                                             N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)

fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh, T_ref_in_evap,
                 T_ref_out_evap, IL, Q_evap, P_evap, T_sat_evap,Le)

dt = 0.01
tspan = (0.0, 1.0)

t = tspan[1]:dt:tspan[2]
len_vec = length(tspan[1]:dt:tspan[2])


ωₐᵢᵣ = zeros(len_vec);
iₐᵢᵣ = zeros(len_vec);
ṁₛₒₗ = zeros(len_vec);
ξₛₒₗ = zeros(len_vec);
iₛₒₗ = zeros(len_vec);

@time simulate!(plateFinCircularTube,fluidThermalData, dt,tspan,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ)

using InteractiveUtils:@code_warntype

@code_warntype simulate!(plateFinCircularTube,fluidThermalData, dt,tspan,ωₐᵢᵣ,iₐᵢᵣ,ṁₛₒₗ,ξₛₒₗ,iₛₒₗ)

@code_warntype NTU(1001.0, 1.0, 0.001, 1e-3, 0.71, 0.002,
            .002, 0.0004, 1.02, 10, 15, 0.002, 14000.1, 1.0, 0.5255,0.025)


using Plots
ωₐᵢᵣ = reverse(ωₐᵢᵣ)
iₐᵢᵣ = reverse(iₐᵢᵣ)

T_air_out_deh = 19.57 + 273.15 # K
RH_out_deh = 0.01 * 49.61 # % 
@show ωₐ_out = HAPropsSI("W", "T", T_air_out_deh, "P", 101325.0 , "R", RH_out_deh)
ωₐᵢᵣ[end]
@show iₐ_out = HAPropsSI("H", "T", T_air_out_deh, "P", 101325.0 , "W", ωₐ_out)
iₐᵢᵣ[end]
plot(t, ωₐᵢᵣ, label = "ωᵢᵣ")
plot(t, iₐᵢᵣ, label = "iᵢᵣ")
plot(t, ṁₛₒₗ, label = "ṁₛₒₗ")
plot(t, ξₛₒₗ, label = "ξₛₒₗ")
plot(t, iₛₒₗ, label = "iₛₒₗ")

ṁₛₒₗ[end] * iₛₒₗ[end]  - ṁₛₒₗ[1] * iₛₒₗ[1]

4 * 0.75 * 0.3 * 7 * 5 * 4.44

(0.005 / 0.2 + 1 / 5) ^ -1

T_sol_calc = (i,ξ) -> calculate_T_sol(i, ξ,IL) - 273.15
T_sol = @. T_sol_calc(iₛₒₗ, ξₛₒₗ)
plot(t, T_sol, label = "Tₛₒₗ")

calculate_T_air = (ω,i) -> HAPropsSI("T", "W", ω , "H", i, "P", 101325.0)
T_air  = @. calculate_T_air(ωₐᵢᵣ,iₐᵢᵣ) - 273.15
plot(t, T_air, label = "Tᵢᵣ")

calculate_T_dp = (ω,i) -> HAPropsSI("Tdp", "W", ω , "H", i, "P", 101325.0)
T_dp  = @. calculate_T_dp(ωₐᵢᵣ,iₐᵢᵣ) - 273.15
plot(t, T_dp, label = "Tdp")

calculate_T_wet = (ω,i) -> HAPropsSI("Twb", "W", ω , "H", i, "P", 101325.0)
T_wet  = @. calculate_T_wet(ωₐᵢᵣ,iₐᵢᵣ) - 273.15
plot(t, T_wet, label = "Twb")

calculate_RH = (ω,i) -> HAPropsSI("R", "W", ω , "H", i, "P", 101325.0)
RH  = @. calculate_RH(ωₐᵢᵣ,iₐᵢᵣ)
plot(t, RH, label = "RH")

i_a_in = HAPropsSI("H", "T", T_air_amb, "P", 101325.0 , "Twb", T_wb_air_amb)
i_a_out = HAPropsSI("H", "T", T_air_out_deh, "P", 101325.0 , "R", RH_out_deh)
ṁₐᵢᵣ_ᵢₙ * (i_a_in - i_a_out) * 2 * N_fin_evap


begin "test evaporator 2 "
    T_air_amb = 25.65 + 273.15  # K
    T_wb_air_amb = 21.2 + 273.15   # K
    m_dot_air_deh = 0.03415  # kg/s
    m_dot_sol_deh = 2.53 / 60 / 1000  # kg/s
    T_sol_in_deh = 19 + 273.15  # K
    X_sol_in_deh = 0.76
    T_ref_in_evap = 9.07 + 273.15  # K
    T_ref_out_evap = 12.86 + 273.15  # K
    P_evap = 3.8 * 100000  # Pa
    T_sat_evap = 6.85 + 273.15  # K
    Q_evap = 1466  # W
    # ========================================
    # T_air_out_deh = 17.54 + 273.15 # K  
    # RH_out_deh = 47.47 # %
    # ========================================
end

begin "test condenser 1 "
    T_air_amb = 25.48 + 273.15  # K
    T_wb_air_amb  = 21.21 + 273.15  # K
    m_dot_air_reg = 0.04725  # kg/s
    m_dot_sol_reg = 2.63 / 60 / 1000  # kg/s
    T_sol_in_reg = 31 + 273.15  # K
    X_sol_in_reg = 0.81
    m_dot_ref = 27.16 / 3600  # kg/s
    T_ref_in_cond = 42.6 + 273.15  # K
    T_ref_out_cond = 37.61 + 273.15  # K
    P_cond = 9.7 * 100000  # Pa
    T_sat_cond = 37.54 + 273.15  # K
    Q_cond = 1192  + 201.7 # W
    # ========================================
    # T_air_out_reg = 32.64 + 273.15  # K
    # RH_out_reg = 51.84  # %
    # ========================================
end

begin "test condenser 2 "
    T_air_amb = 25.65 + 273.15  # K
    T_wb_air_amb = 21.2 + 273.15   # K
    m_dot_air_reg = 0.04731  # kg/s
    m_dot_sol_reg = 2.67 / 60 / 1000  # kg/s
    T_sol_in_reg = 32 + 273.15  # K
    X_sol_in_reg = 0.825
    m_dot_ref = 34.69 / 3600  # kg/s
    T_ref_in_cond = 46.89 + 273.15  # K
    T_ref_out_cond = 40.33 + 273.15  # K
    P_cond = 10.4 * 100000  # Pa
    T_sat_cond = 40.12 + 273.15  # K
    Q_cond = 1466 + 333.5  # W
    # ========================================
    # T_air_out_reg = 34.32 + 273.15  # K 
    # RH_out_reg = 50.61  # %
    # ========================================
end


