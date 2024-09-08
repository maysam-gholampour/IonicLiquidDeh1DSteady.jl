begin 
    δ_fin_evap = 0.00013
    D_tube_outside_evap = 0.0101
    N_tube_per_row_evap = 6
    N_row_evap = 6
    N_fin_evap = 48
    FD_evap = 0.155
    H_evap = 0.132
    FS_evap = 0.00254
    Le = 0.85
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
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh, T_ref_in_evap,
    T_ref_out_evap, IL, Q_evap, P_evap, T_sat_evap,Le)
    # ========================================
    dt = 0.0001
    tspan = (0.0, 1.0)
    t = tspan[1]:dt:tspan[2]
    len_vec = length(tspan[1]:dt:tspan[2])
    # ========================================
end

simulate!(plateFinCircularTube,fluidThermalData, dt,tspan)
