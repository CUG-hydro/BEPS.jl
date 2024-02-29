# LSM of Xuanze Zhang
# 
# Soil moisture is predicted from a 5-layer model (as with soil
# temperature), in which the vertical soil moisture transport is governed
# by infiltration, runoff, gradient diffusion, gravity, and root
# extraction through canopy transpiration.  The net water applied to the
# surface layer is the snowmelt plus precipitation plus the throughfall
# of canopy dew minus surface runoff and evaporation.
# CLM3.5 uses a zero-flow bottom boundary condition.

function UpdateSoilMoisture(p::Soil, kstep::Float64)
  inf, inf_max = 0.0, 0.0
  this_step, total_t, max_Fb = 0.0, 0.0, 0.0

  n = p.n_layer
  @unpack dz, f_water, Ksat, Kb, KK, km, b,
  ψ_sat, ψ, ψb,
  θ_sat, θ, θb, θ_prev = p
  θ_prev .= θ # assign current soil moisture to prev

  # TODO: check this
  @inbounds for i in 1:n+1
    if p.Tsoil_c[i] > 0.0
      f_water[i] = 1.0
    elseif p.Tsoil_c[i] < -1.0
      f_water[i] = 0.1
    else
      f_water[i] = 0.1 + 0.9 * (p.Tsoil_c[i] + 1.0)
    end
  end

  # Max infiltration calculation
  inf_max = f_water[1] * Ksat[1] * (1 + (θ_sat[1] - θ_prev[1]) / dz[1] * ψ_sat[1] * b[1] / θ_sat[1])
  inf = max(f_water[1] * (p.Zp / kstep + p.r_rain_g), 0)
  inf = clamp(inf, 0, inf_max)

  # Ponded water after runoff. This one is related to runoff. LHe.
  p.Zp = (p.Zp / kstep + p.r_rain_g - inf) * kstep * p.r_drainage

  @inbounds while total_t < kstep
    # soil moisture in the boundaries
    for i in 1:n
      if i < n
        θb[i] = (θ[i+1] / dz[i+1] + θ[i] / dz[i]) / (1 / dz[i] + 1 / dz[i+1])
      else
        d1 = max((θ[i] - θb[i-1]) * 2.0 / dz[i], 0)
        θb[i] = min(θ[i] + d1 * dz[i] / 2.0, θ_sat[i])
      end
    end

    # # Kb not used
    # for i in 1:n-1
    #   K = f_water[i] * (θb[i] / θ_sat[i])^(2 * b[i] + 3)
    #   Kb[i] = (Ksat[i] * dz[i] + Ksat[i+1] * dz[i+1]) / (dz[i] + dz[i+1]) * K # bottom
    # end
    # Kb[n] = 0.5 * f_water[n] * cal_K(θ[n], θ_sat[n], Ksat[n], b[n])

    # the unsaturated soil water retention. LHe
    for i in 1:n
      ψ[i] = cal_ψ(θ[i], θ_sat[i], ψ_sat[i], b[i])
      ψb[i] = cal_ψ(θb[i], θ_sat[i], ψ_sat[i], b[i])
    end

    # Hydraulic conductivity: Bonan, Table 8.2, Campbell 1974, K = K_sat*(θ/θ_sat)^(2b+3)
    for i in 1:n
      km[i] = f_water[i] * cal_K(θ[i], θ_sat[i], Ksat[i], b[i])
    end

    # the unsaturated hydraulic conductivity of soil layer @ boundaries
    for i in 1:n-1
      KK[i] = (km[i] * ψ[i] + km[i+1] * ψ[i+1]) / (ψ[i] + ψ[i+1]) * (b[i] + b[i+1]) / (b[i] + b[i+1] + 6)
    end
    # KK[n] = (km[n] * ψ[n] + Kb[n] * ψb[n]) / (ψ[n] + ψb[n]) * b[n] / (b[n] + 3)

    # Fb, flow speed. Dancy's law. LHE.
    # check the r_waterflow further. LHE
    for i in 1:n-1
      Q = KK[i] * (2 * (ψ[i+1] - ψ[i]) / (dz[i] + dz[i+1]) + 1) # z direction
      Q_max = (θ_sat[i+1] - θ[i+1]) * dz[i+1] / kstep + p.Ett[i+1]
      Q = min(Q, Q_max)

      p.r_waterflow[i] = Q
      max_Fb = max(max_Fb, abs(Q))
    end
    # p.r_waterflow[n] = 0

    this_step = guess_step(max_Fb)
    total_t += this_step
    total_t > kstep && (this_step -= (total_t - kstep))

    # from there: kstep is replaced by this_step. LHE
    for i in 1:n
      if i == 1
        θ[i] += (inf - p.r_waterflow[i] - p.Ett[i]) * this_step / dz[i]
      else
        θ[i] += (p.r_waterflow[i-1] - p.r_waterflow[i] - p.Ett[i]) * this_step / dz[i]
      end
      θ[i] = clamp(θ[i], p.theta_vwp[i], θ_sat[i])
    end
  end

  for i in 1:n
    p.ice_ratio[i] *= θ_prev[i] / θ[i]
    p.ice_ratio[i] = min(1.0, p.ice_ratio[i])
  end
end

function guess_step(max_Fb)
  if max_Fb > 1.0e-5
    this_step = 1.0
  elseif max_Fb > 1.0e-6
    this_step = 30.0
  else
    this_step = 360.0
  end
  this_step
end


# Function to compute soil water stress factor
function soil_water_factor_v2(p::Soil)
  # @unpack ψ, ψ_sat, 
  θ = p.θ
  n = p.n_layer

  t1 = -0.02
  t2 = 2.0

  if p.ψ[1] <= 0.000001
    for i in 1:n
      p.ψ[i] = cal_ψ(θ[i], p.θ_sat[i], p.ψ_sat[i], p.b[i])
    end
  end

  for i in 1:n
    # psi_sr in m H2O! This is the old version. LHE.
    p.fpsisr[i] = p.ψ[i] > p.ψ_min ? 1.0 / (1 + ((p.ψ[i] - p.ψ_min) / p.ψ_min)^p.alpha) : 1.0

    p.ft[i] = p.Tsoil_p[i] > 0.0 ? 1.0 - exp(t1 * p.Tsoil_p[i]^t2) : 0
    p.fpsisr[i] *= p.ft[i]
  end

  for i in 1:n
    p.dtt[i] = FW_VERSION == 1 ? p.f_root[i] * p.fpsisr[i] : p.f_root[i]
  end
  dtt_sum = sum(p.dtt)

  if dtt_sum < 0.000001
    p.f_soilwater = 0.1
  else
    fpsisr_sum = 0
    for i in 1:n
      p.dt[i] = p.dtt[i] / dtt_sum
      isnan(p.dt[i]) && println(p.dt[1])

      fpsisr_sum += p.fpsisr[i] * p.dt[i]
    end
    p.f_soilwater = max(0.1, fpsisr_sum)
  end
end


# Function to calculate soil water uptake from a layer
function Soil_Water_Uptake(p::Soil, Trans_o::Float64, Trans_u::Float64, Evap_soil::Float64)
  ρ_w = 1025.0
  Trans = Trans_o + Trans_u

  # for the top layer
  p.Ett[1] = Trans / ρ_w * p.dt[1] + Evap_soil / ρ_w

  # for each layer:
  for i in 2:p.n_layer
    p.Ett[i] = Trans / ρ_w * p.dt[i]
  end
end
