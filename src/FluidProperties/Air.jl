export _μₐ, _kₐ, _cpₐ, _ρₐ, _νₐ, _αₐ

_μₐ(T) = CoolProp.PropsSI("V", "T", T, "P", 101325.0, "Air")
_kₐ(T,ω) = CoolProp.HAPropsSI("K", "T", T, "P", 101325.0, "W", ω)
_cpₐ(T,ω) = CoolProp.HAPropsSI("C", "T", T, "P", 101325.0, "W", ω)
_ρₐ(T,ω) = 1.0 / CoolProp.HAPropsSI("V", "T", T, "P", 101325.0, "W", ω)
_νₐ(T,ω)  = _μₐ(T) / _ρₐ(T,ω) 
_αₐ(T,ω)  = _kₐ(T,ω)  / (_ρₐ(T,ω)  * _cpₐ(T,ω))



