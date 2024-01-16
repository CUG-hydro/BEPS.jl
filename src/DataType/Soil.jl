T = Vector{Float64}

@with_kw mutable struct Soil
  flag        ::Cint    = Cint(0)
  n_layer     ::Cint    = Cint(5)
  step_period ::Cint    = Cint(1)

  Zp          ::Cdouble = Cdouble(0)
  Zsp         ::Cdouble = Cdouble(0)
  r_rain_g    ::Cdouble = Cdouble(0)
  soil_r      ::Cdouble = Cdouble(0)
  r_drainage  ::Cdouble = Cdouble(0)
  r_root_decay::Cdouble = Cdouble(0)
  psi_min     ::Cdouble = Cdouble(0)
  alpha       ::Cdouble = Cdouble(0)
  f_soilwater ::Cdouble = Cdouble(0)

  d_soil      ::T = zeros(10)
  f_root      ::T = zeros(10)
  dt          ::T = zeros(10)
  thermal_cond::T = zeros(10)
  theta_vfc   ::T = zeros(10)
  theta_vwp   ::T = zeros(10)
  fei         ::T = zeros(10)
  Ksat        ::T = zeros(10)
  psi_sat     ::T = zeros(10)
  b           ::T = zeros(10)
  density_soil::T = zeros(10)
  f_org       ::T = zeros(10) # organic matter
  ice_ratio   ::T = zeros(10) # ice ratio
  thetam      ::T = zeros(10) # soil moisture
  thetam_prev ::T = zeros(10) # soil moisture in previous time
  temp_soil_p ::T = zeros(10) # soil temperature in previous time
  temp_soil_c ::T = zeros(10) # soil temperature in current time
  f_ice       ::T = zeros(10) 
  psim        ::T = zeros(10)
  thetab      ::T = zeros(10)
  psib        ::T = zeros(10)
  r_waterflow ::T = zeros(10)
  km          ::T = zeros(10) # hydraulic conductivity
  Kb          ::T = zeros(10)
  KK          ::T = zeros(10) # average conductivity of two soil layers
  Cs          ::T = zeros(10)
  lambda      ::T = zeros(10)
  Ett         ::T = zeros(10)
  G           ::T = zeros(10)

  # temporary variables in soil_water_factor_v2
  ft          ::T = zeros(10)
  dtt         ::T = zeros(10)
  fpsisr      ::T = zeros(10)
end

# ft = zeros(Float64, p.n_layer)
# fpsisr = zeros(Float64, p.n_layer)
# dtt = zeros(Float64, p.n_layer)
