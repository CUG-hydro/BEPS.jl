using Test
using BEPS
import BEPS: kstep


function Base.getindex(x::Union{Soil, Soil_c}, i::Integer)
  names = fieldnames(typeof(x))
  getfield(x, names[i])
end

function is_soil_equal(p_jl, p_c; tol=1e-7, verbose=false)
  names_c = fieldnames(typeof(p_c))
  names_jl = fieldnames(typeof(p_jl))
  ## 变量名可能不同
  
  for i in eachindex(names_c)
    x_jl = getfield(p_jl, i)
    x_c = getfield(p_c, i)

    if verbose
      color = names_c[i] == names_jl[i] ? :black : :red
      printstyled("C and Julia: $(names_c[i]), $(names_jl[i])\n"; color)
    end
    @test maximum(abs.(x_c .- x_jl)) <= tol
  end
end


# Init_Soil_Parameters(par.landcover, par.soil_type, parameter[28], p_soil)
@testset "Init_Soil_Parameters" begin
  p_jl = Soil()
  p_c = Soil_c()

  for stxt = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    Init_Soil_Parameters(1, stxt, 0.1, p_jl)
    clang.Init_Soil_Parameters(1, stxt, 0.1, p_c)
    is_soil_equal(p_jl, p_c; tol=1e-8)
  end
end


@testset "Init_Soil_var_o" begin
  p_jl = Soil()
  p_c = Soil_c()
  var_jl = zeros(41)
  var_c = zeros(41)

  par = (lon=120.5, lat=30.5, landcover=25, clumping=0.85,
    soil_type=8, Tsoil=2.2,
    soilwater=0.4115, snowdepth=0.0)

  rad = 400.0
  temp = 16.0
  hum = 7.0
  pre = 0.0
  wind = 2.0
  meteo = ClimateData(rad, 0.0, temp, hum, pre, wind)

  param = readparam(par.landcover)      # n = 48

  Init_Soil_var_o(p_jl, var_jl, meteo, param, par)
  Init_Soil_var_o(p_c, var_c, meteo, param, par)

  
  funs = [Update_Cs, 
    Update_ice_ratio,
    UpdateSoilThermalConductivity, 
    soil_water_factor_v2]
  
  for fun in funs
    fun(p_jl)
    fun(p_c)
  end

  UpdateSoilMoisture(p_jl, kstep)
  UpdateSoilMoisture(p_c, kstep)

  UpdateHeatFlux(p_jl, 0.0, 0.0, 0.0, 20.0, 3600.0)
  UpdateHeatFlux(p_c, 0.0, 0.0, 0.0, 20.0, 3600.0)

  is_soil_equal(p_jl, p_c; tol=1e-7, verbose=true)
end

# r_root_decay
# f_root
