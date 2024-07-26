using IonicLiquidDeh1DSteady


const N_fin_evap = 48
const N_fin_cond = 75
const FD = 0.205
const g = 9.81
const δₛₒₗ = ∛(3 * ṁₛₒₗ * νₛₒₗ / (ρₛₒₗ * g * FD))
const H = 0.132
const FS = 0.00254
const Reₛₒₗ = Uₛₒₗ_ᵣ * δₛₒₗ / νₛₒₗ




begin "test evaporator 1 "
    T_air_amb = 21.21 + 273.15  # K
    T_wb_air_amb = 25.48 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 2.57 / 60 / 1000  # kg/s
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.76
    T_ref_in_evap = 12.19 + 273.15  # K
    T_ref_out_evap = 17.04 + 273.15  # K
    P_evap = 4.3 * 100000  # Pa
    T_sat_evap = 10.51 + 273.15  # K
    Q_evap = 1192  # W
    # ========================================
    ṁₐ = (m_dot_air_deh / N_fin_evap) * 0.5 # mass flow rate for half of the fin space
    ṁₛₒₗ = (m_dot_sol_deh / N_fin_evap) * 0.5 # mass flow rate for half of the fin space


    # ========================================
    # T_air_out_deh = 19.57 + 273.15 # K
    # RH_out_deh = 49.61 # % 
    # ========================================


end

begin "test evaporator 2 "
    T_air_amb = 21.2 + 273.15   # K
    T_wb_air_amb = 25.65 + 273.15  # K
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
    T_air_amb = 21.21 + 273.15  # K
    T_wb_air_amb = 25.48 + 273.15  # K
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
    T_air_amb = 21.2 + 273.15   # K
    T_wb_air_amb = 25.65 + 273.15  # K
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


