begin
    using Interpolations, StaticArrays
    using NonlinearSolve
    using SimpleNonlinearSolve:ITP
    using CoolProp
    using BoundaryValueDiffEq
    using LinearAlgebra
end

begin "i_vapor_saturation"
    coeffs = @SVector[-7.664856985220677e-5,0.028709375602901342,-6.289410425918317,2207.8158371199142,2.4958805266423156e6]
    function i_vapor_saturation(T)
        T_vec = [T^4, T^3, T^2, T, 1]
        return dot(coeffs, T_vec)
    end
end


begin "Properties Interpolations and Extrapolations"
    const Tⁿᵒᵈᵉˢ = @SVector[x + 273.15 for x in [25.0,35.0,60.0,80.0]]
    const ξⁿᵒᵈᵉˢ_1 = @SVector[x * 0.01 for x in [0.0, 50.0, 70.0, 80.0, 85.0, 90.0]]
    const ξⁿᵒᵈᵉˢ_2 = @SVector[x * 0.01 for x in [0.0, 50.0, 70.0, 80.0, 85.0, 90.0, 95.0]]
    const nodes_P_ν = (Tⁿᵒᵈᵉˢ, ξⁿᵒᵈᵉˢ_1)
    const nodes = (Tⁿᵒᵈᵉˢ, ξⁿᵒᵈᵉˢ_2)
    const P_v_data = @SMatrix[
        3200.0   2600.0   1700.0   1000.0    600.0   250.0
        5600.0   4500.0   3000.0   1800.0   1200.0   500.0
        20000.0  16300.0  11200.0   7100.0   4900.0  2600.0
        47500.0  39300.0  27900.0  18200.0  12800.0  7100.0
        ]
    const ρ_data = @SMatrix[  
        997.0  1092.0  1128.0  1143.0  1147.0  1151.0  1154.0;
        994.0  1086.0  1121.0  1136.0  1140.0  1144.0  1147.0;
        983.0  1068.0  1103.0  1117.0  1121.0  1125.0  1128.0;
        972.0  1054.0  1088.0  1102.0  1106.0  1110.0  1113.0
        ]
    const μ_data = @SMatrix[
        0.001   0.01   0.02   0.032  0.044  0.062  0.105
        0.0007  0.006  0.013  0.021  0.027  0.039  0.063
        0.0005  0.004  0.006  0.009  0.011  0.015  0.023
        0.0003  0.003  0.004  0.005  0.006  0.008  0.012
        ]
    const cp_data = @SMatrix[
        4200.0  2900.0  2400.0  2200.0  2100.0  1900.0  1700.0
        4200.0  2900.0  2400.0  2200.0  2100.0  1900.0  1700.0
        4200.0  3000.0  2500.0  2300.0  2200.0  2000.0  1900.0
        4200.0  3100.0  2600.0  2400.0  2300.0  2100.0  2000.0
        ]
    const 𝑘_data = @SMatrix[
        0.607  0.377  0.301  0.274  0.266  0.257  0.249;
        0.62   0.381  0.304  0.276  0.268  0.258  0.251;
        0.65   0.393  0.313  0.281  0.273  0.263  0.256;
        0.67   0.403  0.32   0.285  0.277  0.267  0.259
        ]
    const i_data = @SMatrix[
        0.0  -58000.0  -75000.0  -74000.0  -68000.0  -55000.0  -34000.0
        0.0  -57000.0  -72000.0  -72000.0  -67000.0  -54000.0  -33000.0
        0.0  -52000.0  -67000.0  -67000.0  -62000.0  -51000.0  -31000.0
        0.0  -48000.0  -62000.0  -64000.0  -59000.0  -49000.0  -30000.0
        ]
    # ================== Interpolation and Extrapolation P_ν ==================
    P_ν_interplated = interpolate(nodes_P_ν, P_v_data, Gridded(Linear()))
    P_ν_extrapolated = extrapolate(P_ν_interplated, Line())
    function Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ)
        return max(P_ν_extrapolated(T, ξ),0.0001)
    end
    # ================== Interpolation and Extrapolation ρ ====================
    ρ_interplated = interpolate(nodes, ρ_data, Gridded(Linear()))
    ρ_extrapolated = extrapolate(ρ_interplated, Line())
    function ρₛₒₗ(T, ξ)
        return ρ_extrapolated(T, ξ)
    end
    # ================== Interpolation and Extrapolation μ ====================
    μ_interplated = interpolate(nodes, μ_data, Gridded(Linear()))
    μ_extrapolated = extrapolate(μ_interplated, Line())
    function μₛₒₗ(T, ξ)
        return μ_extrapolated(T, ξ)
    end
    # ================== Interpolation and Extrapolation cp ====================
    cp_interpolated = interpolate(nodes, cp_data, Gridded(Linear()))
    cp_extrapolated = extrapolate(cp_interpolated, Line())
    function cpₛₒₗ(T, ξ)
        return cp_extrapolated(T, ξ)
    end
    # ================== Interpolation and Extrapolation 𝑘 ====================
    𝑘_interpolated = interpolate(nodes, 𝑘_data, Gridded(Linear()))
    𝑘_extrapolated = extrapolate(𝑘_interpolated, Line())
    function 𝑘ₛₒₗ(T, ξ)
        return 𝑘_extrapolated(T, ξ)
    end
    # ================== Interpolation and Extrapolation 𝑖 ====================
    i_interpolated = interpolate(nodes, i_data, Gridded(Linear()))
    i_extrapolated = extrapolate(i_interpolated, Line())
    function iₛₒₗ(T, ξ)
        return i_extrapolated(T, ξ)
    end
    # ================== Find T given i_sol and ξ ====================
    # Function to find the root, given i_sol and ξ
    function calculate_T(iᵛₛₒₗ, ξ; T_lower=0.0 + 273.15, T_upper=95.0 + 273.15)
        f(T, p) = iₛₒₗ(T, p[2]) - p[1]
        p = SA[iᵛₛₒₗ,ξ]
        T_span = @SVector[T_lower , T_upper]
        prob = IntervalNonlinearProblem(f, T_span, p)
        result = solve(prob, ITP())
        return calculate_T_barrier(result)
    end
    function calculate_T_barrier(result)
        return result.u
    end
end


begin "test Properties"
    # ================== Test ====================
    # @code_warntype P_ν_sol(25.0, 0.1); # 3.2
    @show Pᵥₐₚₒᵣ_ₛₒₗ(25.0 + 273.15, 0.1 * 0.01); # 3.2
    @show ρₛₒₗ(25.0 + 273.15, 0.1 * 0.01); # 997.0
    @show μₛₒₗ(25.0 + 273.15, 0.1 * 0.01); # 1.0
    @show cpₛₒₗ(25.0 + 273.15, 0.1 * 0.01); # 4.2
    @show 𝑘ₛₒₗ(25.0 + 273.15, 0.1 * 0.01); # 0.607
    @show iₛₒₗ(25.0 + 273.15, 2 * 0.01); # 0.0

    @show i = iₛₒₗ(25.0 + 273.15, 74.0 * 0.01) 

    # @code_warntype find_T(𝑖, 74.0); # 25.0
    @code_warntype calculate_T(i, 74.0 * 0.01) - 273.15
end


function ionic_liquid_coil_ode!(du,u, p, t)
    # ωₐᵢᵣ, iₐᵢᵣ, ṁₛₒₗ,ξₛₒₗ, iₛₒₗ = u
    # Pₐᵢᵣ, H, Le, ∂Qᵣ , ṁₐᵢᵣ, NTUᴰₐᵢᵣ, σ = p
    Tₛₒₗ = calculate_T(u[5],u[4])
    Pᵥₐₚₒᵣ_ₛₒₗᵛ = Pᵥₐₚₒᵣ_ₛₒₗ(Tₛₒₗ,u[4])
    ωₑ = 0.622 * Pᵥₐₚₒᵣ_ₛₒₗᵛ / (p[1] - Pᵥₐₚₒᵣ_ₛₒₗᵛ)
    iₑ = 1.01 * (Tₛₒₗ - 273.15) + ωₑ * (2500 + 1.04 * (Tₛₒₗ - 273.15))
    iₑ *= 1000
    # iᵥₐₚₒᵣ_ₜₛ = PropsSI("H", "T", Tₛₒₗ + 273.15, "Q", 1, "water")
    iᵥₐₚₒᵣ_ₜₛ = i_vapor_saturation(Tₛₒₗ)

    du[1] = p[7] * (p[6] / p[2]) * (u[1] - ωₑ)
    du[2] = p[7] * (p[6] * p[3] / p[2]) * ((u[2] - iₑ) + (iᵥₐₚₒᵣ_ₜₛ * (1 / p[3] - 1) * (u[1] - ωₑ)))
    du[3] = p[7] * p[5] * du[1]
    du[4] = (-u[4] / u[3]) * du[3]
    du[5] = (1 / u[3]) * (p[7] * p[5] * du[2] - u[5] * du[3] - p[4])
    nothing
end
# ω_end =  0.01782, i_end = 75725.5949 air
# ṁs0= 0.0006004891705935622, ξ0= 0.7, i0 = -74600.0
function bca!(res_a, u_a, p)
    res_a[1] = u_a[3] - 0.0006004891705935622
    res_a[2] = u_a[4] - 0.74
    res_a[3] = u_a[5] + 74600.0
    nothing
end
function bcb!(res_b, u_b, p)
    res_b[1] = u_b[1] - 0.01782
    res_b[2] = u_b[2] - 75725.5949
    nothing
end

p = @SVector[101325.0, 0.135, 0.8, 250, 0.0007506114632419527, 3.125, 1.0];
u0 = [0.009, 65000.5949271956, 0.0006004891705935622 + 0.00001, 0.65, -74600.0];
tspan = (0.0, 0.15)

bvp_fun = BVPFunction(
    ionic_liquid_coil_ode!, (bca!, bcb!),
    bcresid_prototype = (zeros(3), zeros(2)), twopoint = Val(true))


prob = BVProblem(bvp_fun,
    u0,
    tspan,
    p)

@time sol = solve(prob, MIRK4(), dt = 0.001);

prob_new = remake(prob, p=@SVector[101325.0, 0.135, 0.8, 200, 0.0007506114632419527, 3.125, 1.0];)
@time sol_new = solve(prob_new, MIRK4(), dt = 0.001);

size(sol.u[1])
x1 = zeros(length(sol.u))
x2 = zeros(length(sol.u))
x3 = zeros(length(sol.u))
x4 = zeros(length(sol.u))
x5 = zeros(length(sol.u))
for i in 1:length(sol.u)
    x1[i] = sol.u[i][1]
    x2[i] = sol.u[i][2]
    x3[i] = sol.u[i][3]
    x4[i] = sol.u[i][4]
    x5[i] = sol.u[i][5]
end

using Plots
x1 = reverse(x1)
x2 = reverse(x2)
plot(sol.t, x1, label = "ωᵢᵣ")
plot(sol.t, x2, label = "iᵢᵣ")
plot(sol.t, x3, label = "ṁₛₒₗ")
plot(sol.t, x4, label = "ξₛₒₗ")
plot(sol.t, x5, label = "iₛₒₗ")

T_sol = @. calculate_T(x5, x4) - 273.15
plot(sol.t, T_sol, label = "Tₛₒₗ")

calculate_T_air = (ω,i) -> HAPropsSI("T", "W", ω , "H", i, "P", 101325.0)
T_air  = @. calculate_T_air(x1,x2) - 273.15
plot(sol.t, T_air, label = "Tᵢᵣ")

calculate_T_dp = (ω,i) -> HAPropsSI("Tdp", "W", ω , "H", i, "P", 101325.0)
T_dp  = @. calculate_T_dp(x1,x2) - 273.15
plot(sol.t, T_dp, label = "Tdp")

calculate_T_wet = (ω,i) -> HAPropsSI("Twb", "W", ω , "H", i, "P", 101325.0)
T_wet  = @. calculate_T_wet(x1,x2) - 273.15
plot(sol.t, T_wet, label = "Twb")

calculate_RH = (ω,i) -> HAPropsSI("R", "W", ω , "H", i, "P", 101325.0)
RH  = @. calculate_RH(x1,x2)
plot(sol.t, RH, label = "RH")