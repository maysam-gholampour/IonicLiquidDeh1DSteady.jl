

begin "Properties Interpolations and Extrapolations"
    
    const Tⁿᵒᵈᵉˢ = @SVector[x + 273.15 for x in [25.0,35.0,60.0,80.0]]
    const ξⁿᵒᵈᵉˢ_2 = @SVector[x * 0.01 for x in [0.0, 50.0, 70.0, 80.0, 85.0, 90.0, 95.0]]
    const nodes = (Tⁿᵒᵈᵉˢ, ξⁿᵒᵈᵉˢ_2)
    const Tⁿᵒᵈᵉˢ_s = @SVector[x + 273.15 for x in [25.0,35.0,60.0]]
    const nodes_s = (Tⁿᵒᵈᵉˢ_s, ξⁿᵒᵈᵉˢ_2)
    const Δh_data = @SMatrix[
        0.0  -58000.0  -75000.0  -74000.0  -68000.0  -55000.0  -34000.0
        0.0  -57000.0  -72000.0  -72000.0  -67000.0  -54000.0  -33000.0
        0.0  -52000.0  -67000.0  -67000.0  -62000.0  -51000.0  -31000.0
        0.0  -48000.0  -62000.0  -64000.0  -59000.0  -49000.0  -30000.0
        ]
    const σ_data = @SMatrix[
        0.0719  0.0441  0.0399  0.0388  0.0376  0.0367  0.036
        0.0704  0.043   0.0391  0.0382  0.0371  0.0364  0.0359
        0.0662  0.0403  0.0371  0.0366  0.0357  0.0351  0.0347
        ]
    # ================== Interpolation and Extrapolation P_ν ==================
    const a0_p = 12.10 
    const a1_p = -28.01 
    const a2_p = 50.34 
    const a3_p = -24.63
    const b0_p = 1212.67 
    const b1_p = 772.37 
    const b2_p = 614.59 
    const b3_p = 493.33

    # @inline function _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ,::CreCOPlus5100)
    #     A = a0_p + a1_p * ξ + a2_p * ξ^2 + a3_p * ξ^3
    #     B = b0_p + b1_p * ξ + b2_p * ξ^2 + b3_p * ξ^3
    #     return 10^(A - B / T) * 100.0
    # end #original function

    @inline function _Pᵥₐₚₒᵣ_ₛₒₗ(T, ξ,::CreCOPlus5100)
        A = a0_p + a1_p * ξ + a2_p * ξ^2 + a3_p * ξ^3
        B = b0_p + b1_p * ξ + b2_p * ξ^2 + b3_p * ξ^3
        return 10^(A - B / T) * 100.0 * 1.03
    end
    # ================== Interpolation and Extrapolation ρ ====================
    function _ρₛₒₗ(T, ξ,::CreCOPlus5100)
        a0 = 804.28 + 1.585 * T - 0.0031 * T^2
        a1 = 1036.04 - 4.42 * T + 0.0057 * T^2
        a2 = -403.62 + 1.745 * T - 0.0021  * T^2
        return a0 + a1 * ξ + a2 * ξ^2
    end
    # ================== Interpolation and Extrapolation μ ====================
    function _μₛₒₗ(T, ξ,::CreCOPlus5100)
        return (10 ^ (-1.315 + 6.1044 * ξ + 1.0944 * ξ^2 + 342.2 / T - 0.0163 * ξ * T)) * 1e-3
    end
    # ================== Interpolation and Extrapolation cp ====================
    
    function _cpₛₒₗ(T, ξ,::CreCOPlus5100)
        return ((0.00476 * T - 4.01)* ξ + 4.21) * 1e3 
    end
    # ================== Interpolation and Extrapolation 𝑘 ====================

    function _𝑘ₛₒₗ(T, ξ,::CreCOPlus5100)
        return (-0.00103* T - 0.077) * ξ + 0.00109 * T + 0.267
    end
    # ================== Interpolation and Extrapolation 𝑖 ====================
    # @inline function _Δh(T, ξ,::CreCOPlus5100)
    #     p1 = 83.57264585405902
    #     p2 = 55.23823946168705
    #     p3 = 53.55061728738012
    #     p4 = 83.42362593623382
    #     p5 = -819.8792814164431
    #     p6 = 591.1712249511673
    #     T_c_H2O = 647.226
    #     θ = T / T_c_H2O
    #     Δhₜ = p1 + p2 * θ + p3 * θ ^ 2
    #     Δhₓ = p4 + p5 * ξ  + p6 * ξ ^ 2 
    #     return (Δhₜ + Δhₓ + 1.0) * 1e3
    # end

    @inline function _Δh(T, ξ,::CreCOPlus5100)
        Δh_interpolated = interpolate(nodes, Δh_data, Gridded(Linear()))
        Δh_extrapolated = extrapolate(Δh_interpolated, Line())
        return Δh_extrapolated(T, ξ)
    end

    @inline function _iₛₒₗ(T, ξ,::CreCOPlus5100)
        Δh = _Δh(T, ξ,CreCOPlus5100())
        i = _cpₛₒₗ(T, ξ,CreCOPlus5100()) * (T - 273.15) + Δh
        return i
    end
    # ================== Interpolation and Extrapolation σ ====================
    
    function _σₛₒₗ(T, ξ,::CreCOPlus5100)
        σ_interpolated = interpolate(nodes_s, σ_data, Gridded(Linear()))
        σ_extrapolated = extrapolate(σ_interpolated, Line())
        return σ_extrapolated(T, ξ)
    end
    # ================== Find T given i_sol and ξ ====================
    # Function to find the root, given i_sol and ξ
    @inline function calculate_T_sol(iᵛₛₒₗ, ξ,::CreCOPlus5100 ;T_lower=-150.0 + 273.15, T_upper=95.0 + 273.15) 
        f(T, p)= _iₛₒₗ(T, p[2],CreCOPlus5100()) - p[1]
        p = @SVector[iᵛₛₒₗ,ξ]
        T_span = @SVector[T_lower , T_upper]
        prob = IntervalNonlinearProblem{false}(f, T_span, p)
        result = solve(prob, ITP())
        return result.u
    end
end

