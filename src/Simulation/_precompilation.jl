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
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 0.2 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.8
    Q_evap = 1192.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.05
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation CreCOPlus5100 is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan)
@info "Precompilation CreCOPlus5100 completed"

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
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 1.0 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.8
    Q_evap = 1192.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.001
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation CreCOPlus5100 nested is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK4())
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK5())
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK6())
@info "Precompilation CreCOPlus5100 nested completed"

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
    IL = LiCl()
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 1.0 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.4
    Q_evap = 0.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.005
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation LiCl is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan)
@info "Precompilation LiCl completed"

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
    IL = CaCl2()
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 0.2 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.4
    Q_evap = 0.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.05
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation CaCl2 is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan)
@info "Precompilation CaCl2 completed"

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
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 1.0 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.8
    Q_evap = -1192.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.001
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation CreCOPlus5100 nested is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK4())
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK5())
simulate(plateFinCircularTube,fluidThermalData, dt,tspan,solver = MIRK6())
@info "Precompilation CreCOPlus5100 nested completed"

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
    T_air_amb = 30.0 + 273.15 # K  #NOTE: 1.5 is the Blower temperature rise
    T_wb_air_amb = 27.7 + 273.15  # K
    m_dot_air_deh = 0.03584  # kg/s
    m_dot_sol_deh = 0.2 * m_dot_air_deh  # kg/s  #FIXME: 1100 is density
    T_sol_in_deh = 21 + 273.15  # K
    X_sol_in_deh = 0.8
    Q_evap = -1192.0
    # ========================================
    plateFinCircularTube = PlateFinCircularTube(δ_fin_evap, D_tube_outside_evap, N_tube_per_row_evap, N_row_evap,
    N_fin_evap, FD_evap, H_evap, FS_evap, σ_evap)
    fluidThermalData = FluidThermalData(T_air_amb, T_wb_air_amb, m_dot_air_deh, m_dot_sol_deh, T_sol_in_deh, X_sol_in_deh,IL, Q_evap,Le)
    # ========================================
    dt = 0.05
    tspan = (0.0, 1.0)
    # ========================================
end
@info "Precompilation CreCOPlus5100 is running..."
simulate(plateFinCircularTube,fluidThermalData, dt,tspan)
@info "Precompilation CreCOPlus5100 completed"
