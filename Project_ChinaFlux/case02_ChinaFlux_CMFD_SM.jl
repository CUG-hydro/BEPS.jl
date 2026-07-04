## 1. 读取驱动数据
## CMFD V4 + SM-aware calibration: keeps the original case02_ChinaFlux_CMFD.jl untouched.
using BEPS, RTableTools, DataFrames, Dates, ModelParams, Ipaper, JLD2, Statistics
FT = Float64

cmfd_time(t::DateTime) = t
cmfd_time(t) = DateTime(first(String(t), 19), dateformat"yyyy-mm-ddTHH:MM:SS")
cmfd_local_time(t) = cmfd_time(t) + Hour(8)

# 元数据为手工维护：含 "NA" 等会把整列读成 String，统一转数值（不可解析→missing）
as_num(x) = x isa Number ? Float64(x) : (x isa AbstractString ? something(tryparse(Float64, x), missing) : missing)

env_bool(name, default=false) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "y")
OVERWRITE = env_bool("BEPS_OVERWRITE")
MAXN = parse(Int, get(ENV, "BEPS_MAXN", "1000"))
SM_GOAL = Symbol(get(ENV, "BEPS_SM_GOAL", "KGE"))

function valid_mean(x; fallback=NaN)
  vals = Float64[]
  for v in skipmissing(x)
    y = Float64(v)
    isfinite(y) && y != -999.0 && push!(vals, y)
  end
  isempty(vals) ? fallback : mean(vals)
end

base_paths() = [
  [:r_drainage],
  [:veg, :Ω],
  [:veg, :g1_w],
  [:veg, :g0_w],
  [:veg, :VCmax25],
]

function load_cmfd_initial!(model, SITE)
  f = joinpath(BASE_OUTDIR, "BEPS_$(SITE).jld2")
  isfile(f) || return false
  theta = load(f, "theta_opt")
  length(theta) == length(base_paths()) || return false
  BEPS.update!(model, base_paths(), theta)
  true
end

function sm_layer_indices(model, depths_SM)
  centers = cumsum(model.dz) .- model.dz ./ 2
  sort!(unique([argmin(abs.(centers .- d)) for d in depths_SM]))
end

function sm_paths(model, depths_SM)
  layers = sm_layer_indices(model, depths_SM)
  theta_paths = [[:hydraulic, :profile, :θ_sat, i] for i in layers]
  [theta_paths; [
    [:hydraulic, :profile, :ψ_sat, 1],
    [:hydraulic, :profile, :b, 1],
    [:hydraulic, :kv, :kv, 1],
    [:hydraulic, :kv, :f, 1],
  ]]
end

function fit_kvexp(dz, K_sat)
  z = 100 .* (cumsum(dz) .- dz ./ 2)
  y = log.(max.(K_sat, eps(eltype(K_sat))))
  zbar, ybar = mean(z), mean(y)
  slope = sum((z .- zbar) .* (y .- ybar)) / sum((z .- zbar) .^ 2)
  f = clamp(-slope, 0.0, 0.1)
  kv = clamp(exp(ybar + f * zbar), 0.002, 100.0)
  kv, f
end

function with_kvexp(model::ParamBEPS{FT,N}) where {FT,N}
  kv0, f0 = fit_kvexp(model.dz, model.hydraulic.K_sat)
  kv = KvExpLayers{FT,N}(; kv=fill(FT(kv0), N), f=fill(FT(f0), N))
  hydraulic = HydraulicProfile{FT,N}(deepcopy(model.hydraulic.profile), kv, FT.(100 .* model.dz))
  hydraulic.profile.K_sat .= kv0 .* exp.(-f0 .* (cumsum(hydraulic.dz_cm) .- hydraulic.dz_cm ./ 2))
  ParamBEPS{FT,N}(hydraulic, deepcopy(model.thermal);
    dz=copy(model.dz), r_drainage=model.r_drainage, ψ_min=model.ψ_min,
    alpha=model.alpha, veg=deepcopy(model.veg))
end

function sync_sm_profile!(model)
  n = model.N
  model.hydraulic.profile.ψ_sat[2:n] .= model.hydraulic.profile.ψ_sat[1]
  model.hydraulic.profile.b[2:n] .= model.hydraulic.profile.b[1]
  model.hydraulic.kv.kv[2:n] .= model.hydraulic.kv.kv[1]
  model.hydraulic.kv.f[2:n] .= model.hydraulic.kv.f[1]
  z = cumsum(model.hydraulic.dz_cm) .- model.hydraulic.dz_cm ./ 2
  model.hydraulic.profile.K_sat .= model.hydraulic.kv.kv[1] .* exp.(-model.hydraulic.kv.f[1] .* z)
  nothing
end

function eval_sm(theta::Vector{FT}, model::ParamBEPS{FT}, forcing::MetSeries{FT},
  lai::Vector{FT}, dates_UTC::Vector{DateTime}; paths, lon::FT, lat::FT,
  depths_SM::Vector{FT}, depths_TS::Vector{FT}, FluxDay::DataFrame,
  SolveSM_fn=SolveSM_BEPS, ignored...) where {FT<:AbstractFloat}

  params = parameters(model; paths)
  theta_prev = params.value
  BEPS.update!(model, paths, theta)
  sync_sm_profile!(model)
  state = InitState0(model, forcing)
  df_fluxes, _, states, _ = simulate(forcing, lai, dates_UTC;
    ps=model, state, lon, lat, SolveSM_fn)
  gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_UTC .+ Hour(8), FluxDay;
    depths_SM, depths_TS)
  BEPS.update!(model, paths, theta_prev)
  sync_sm_profile!(model)
  gof, data_sim, data_obs
end

function loss_sm(theta::Vector{FT}, model::ParamBEPS{FT},
  forcing::MetSeries{FT}, lai::Vector{FT}, dates_UTC::Vector{DateTime};
  paths, sm_goal=SM_GOAL, kw_loss...) where {FT<:AbstractFloat}

  gof, _, _ = eval_sm(theta, model, forcing, lai, dates_UTC; paths, kw_loss...)
  hasproperty(gof.SM, sm_goal) || return FT(999)
  -valid_mean(gof.SM[!, sm_goal]; fallback=-999.0)
end

function optim_sm_only(model::ParamBEPS{FT}, forcing::MetSeries{FT}, lai::Vector{FT},
  dates_UTC::Vector{DateTime}; paths, maxn=200, kw_loss...) where {FT<:AbstractFloat}

  params = parameters(model; paths)
  lb = map(x -> FT(x[1]), params.bound)
  ub = map(x -> FT(x[2]), params.bound)
  theta, _, _ = sceua(theta -> loss_sm(theta, model, forcing, lai, dates_UTC; paths, kw_loss...),
    params.value, lb, ub; maxn, verbose=true, parallel=true)
  theta
end

##
const BASE_OUTDIR = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE_CMFD_1h_V4"
outdir = "Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE_CMFD_1h_V4_SM"
indir = "/mnt/z/GitHub/jl-pkgs/ChinaFlux2026"
st_full = fread("$indir/data/Metadata/ChinaFlux_Metadata.csv")  # 39 站元数据，含 31 站全部

f = "/mnt/z/China/CMFD_V2.0/OUTPUT/ChinaFlux/ChinaFlux_sp38_final_forcing_1h.csv"
@time FORCING = fread(f)
replace_missing!(FORCING)
rename_existing!(d, pairs) = rename!(d, filter(p -> first(p) in propertynames(d), pairs))
rename_existing!(FORCING, [
  :datetime => :time,
  :Temp => :Tair,
  :RHum => :RH,
  :Wind => :Uz,
  :SRad => :Rs,
  :LRad => :Rln_in,
  :Prec => :Prcp,
])

flux_files = filter(f -> endswith(f, "_Daily_FluxALL_v20260615.csv"),
  readdir("$indir/data/BEPS/Daily_FluxALL"))
SITES_obs = replace.(flux_files, "_Daily_FluxALL_v20260615.csv" => "")
SITES = sort(intersect(unique(FORCING.site), SITES_obs))

miss_f = setdiff(unique(FORCING.site), SITES_obs)   # forcing 有、观测无
miss_o = setdiff(SITES_obs, unique(FORCING.site))   # 观测有、forcing 无
isempty(miss_f) || @warn "CMFD forcing without daily obs" miss_f
isempty(miss_o) || @warn "Daily obs without CMFD forcing" miss_o

HAS_SITE_FILTER = haskey(ENV, "BEPS_SITES")
if HAS_SITE_FILTER
  SITES = intersect(SITES, strip.(split(ENV["BEPS_SITES"], ",")))
  isempty(SITES) && error("BEPS_SITES 无交集")
end

# V3: 仅重跑冠层高度 > 10m 的站点（之前的 safe_wind_ref 处理有误）
h_overs = Dict(r.site => as_num(r.z_overstory) for r in eachrow(st_full) if !ismissing(as_num(r.z_overstory)))

# SITES = intersect(SITES, [s for (s, h) in h_overs if h > 10])
# @info "Tall canopy sites (h > 10m) to rerun" SITES

# V4: 仅重跑冠层高度 < 10m 的站点，不使用 Rln_in（模型用气温估算长波辐射）
SITES = intersect(SITES, [s for (s, h) in h_overs if h < 10])
@info "Short canopy sites (h < 10m, no Rln_in) to rerun" SITES

# 鲁棒列处理：仅重命名存在的列；保证列为 Float64（缺测→NaN），整列不存在→全 NaN
ascol!(d, col) = d[!, col] = (col in propertynames(d)) ?
                             Float64.(coalesce.(d[!, col], NaN)) : fill(NaN, nrow(d))

##
function LoadData(SITE)
  f = "$indir/data/BEPS/Daily_FluxALL/$(SITE)_Daily_FluxALL_v20260615.csv" |> path_mnt

  FluxALL = fread(f)
  replace_missing!(FluxALL)
  # 部分站缺碳通量列（GPP/Hs），rename 仅作用于存在的列
  rename_existing!(FluxALL, [:LAI_glass_G005 => :lai, :GPP => :GPP_obs, :ET => :ET_obs, :Hs => :Hs_obs])
  foreach(c -> ascol!(FluxALL, c), [:lai, :GPP_obs, :ET_obs, :Hs_obs])  # 缺列填 NaN，GOF 自然返回 NaN
  normalize_flux_obs!(FluxALL)
  # 观测已为标准日尺度单位（GPP: gC m⁻² d⁻¹, ET: mm d⁻¹, Hs: W m⁻²），与 agg_daily 输出一致，无需换算
  # （单位见 data/BEPS/BEPS_Forcing_China_FluxALL.md §1.2）
  (; lai) = FluxALL
  ntime2 = length(lai) * 24

  d_forcing = FORCING[FORCING.site.==SITE, :]
  sort!(d_forcing, :time)

  # 驱动可能比观测更早开始/更长（数据本身如此）：对齐到观测首日，再裁剪到观测长度
  i_beg = findfirst(t -> Date(cmfd_local_time(t)) == FluxALL.date[1], d_forcing.time)
  isnothing(i_beg) && error("forcing 不含观测首日 $(FluxALL.date[1])")
  d_forcing = d_forcing[i_beg:min(i_beg + ntime2 - 1, end), :]

  clean_stats = sanitize_forcing!(d_forcing)
  @info "Forcing quality control" clean_stats
  (; Tair, RH, Uz, Rs, Rln_in, Prcp) = d_forcing

  # V4: 不使用 CMFD 长波辐射，令模型用 cal_Rln(ϵ_air, Tair) 估算（netRadiation.jl:190）
  Rln_in = fill(NaN, length(Tair))
  ntime = length(Tair)
  forcing = MetSeries(; ntime, Rs, Rln_in, Tair, RH, Uz, Prcp)
  dates_local = cmfd_local_time.(d_forcing.time)

  dates_local, forcing, lai, FluxALL
end


function RunModel(SITE; maxn=1000, outdir="Project_ChinaFlux/OUTPUT/ALL/Bonan/NSE_CMFD_1h_V4_SM", overwrite=false,
  goal=:NSE, goal_multiplier=-1, sm_goal=SM_GOAL, SolveSM_fn=SolveSM_Bonan)

  mkpath(outdir)
  fout = "$outdir/BEPS_$(SITE).jld2"
  (isfile(fout) && !overwrite) && return

  printstyled("[site]: $SITE\n", color=:blue, bold=true, underline=true)
  t_beg = time()  # 记录单站运行时长（主要由 optim 决定），存入结果

  # 率定所需数据
  dates_local, forcing, lai, FluxALL = LoadData(SITE)
  dates_UTC = dates_local .- Hour(8) # [local] -> [UTC]

  ## 2. 初始化模型参数和状态变量
  i_st = findfirst(st_full.site .== SITE)
  isnothing(i_st) && error("missing metadata: $SITE")
  st = st_full[i_st, :]
  (; lon, lat, VegType, SoilType) = st
  VegType == "GRA" && (VegType = "CRO")
  z_Uz, z_overstory = as_num(st.z_Uz), as_num(st.z_overstory)  # 元数据手工维护，统一转数值
  SoilType = ismissing(SoilType) ? "loam" : SoilType    # 元数据缺失兜底

  # 元数据 z_SM/z_TS 可能缺失；缺失→不评价土壤
  parse_depths(s) = (ismissing(s) || isempty(strip(String(s)))) ?
                    Float64[] : map(x -> parse(Int, strip(x)), split(String(s), ",")) ./ 100
  depths_SM = parse_depths(st.z_SM)
  depths_TS = parse_depths(st.z_TS)

  # 深度过滤：各站日文件土壤列高度异构、命名可能不规范（如固城 z_SM=4 但列名为 SM_4cm_N）。
  # 仅保留标准列 SM_{d}cm / TS_{d}cm 存在「且含有限观测」的深度；无匹配则该站只评价 GPP/ET/Hs。
  function has_obs(prefix, d)
    c = Symbol("$(prefix)_$(Int(round(d * 100)))cm")
    c in propertynames(FluxALL) && any(isfinite, coalesce.(FluxALL[!, c], NaN))
  end
  depths_SM = filter(d -> has_obs("SM", d), depths_SM)
  depths_TS = filter(d -> has_obs("TS", d), depths_TS)

  model = ParamBEPS(VegType, SoilType) |> with_kvexp
  ismissing(z_overstory) || (model.veg.z_canopy_o = z_overstory)
  # CMFD 风速实测高度为 10 m。src 中 safe_wind_ref 自动处理高冠层的参考高度抬升，
  # 无需调用方手动外推风速。
  model.veg.z_wind = 10.0

  load_cmfd_initial!(model, SITE) || @warn "CMFD 初值缺失，使用默认参数作为 SM 初值" SITE
  isempty(depths_SM) && return

  state = InitState0(model, forcing)
  @time df_fluxes, df_ET, states, caches = simulate(forcing, lai, dates_UTC;
    ps=model, state, lon, lat, SolveSM_fn)
  @time gof, data_sim, data_obs = BEPS_GOF(df_fluxes, states, dates_local, FluxALL;
    depths_SM, depths_TS)

  ## 3. 仅优化 SM 关键水力参数：观测土层 θ_sat + kv
  paths = sm_paths(model, depths_SM)
  isempty(paths) && return
  opts = DataFrame(parameters(model; paths))

  kw_loss = (; lon, lat, depths_SM, depths_TS, FluxDay=FluxALL,
    goal, goal_multiplier, sm_goal, SolveSM_fn)

  @time theta_opt = optim_sm_only(model, forcing, lai, dates_UTC; paths, maxn, kw_loss...)
  opts.theta_opt = theta_opt
  gof_opt, data_sim, data_obs = eval_sm(theta_opt, model, forcing, lai, dates_UTC; paths, kw_loss...)

  runtime = round(time() - t_beg, digits=1)  # [秒]
  jldsave(fout; gof_opt, gof, theta_opt, opts, data_sim, data_obs, runtime)
  gof_opt
end


function Process(SITES)
  errors = Tuple{String,String}[]
  for i in eachindex(SITES)
    !isCurrentWorker(i) && continue

    SITE = SITES[i]
    try
      RunModel(SITE; maxn=MAXN, outdir, goal=:NSE, goal_multiplier=-1,
        SolveSM_fn=SolveSM_Bonan, overwrite=OVERWRITE)
    catch ex
      msg = sprint(showerror, ex)
      @error "Error processing site $SITE: $msg"
      push!(errors, (SITE, msg))
    end
  end

  if !isempty(errors)
    @warn "以下站点运行失败" errors
    fwrite(DataFrame(site=first.(errors), error=last.(errors)), "$outdir/_errors.csv")
  end
end

SITES_poor = [
  "CRO_冬小麦夏玉米_固城",
  "CRO_水稻_盘锦",
  "CRO_水稻_长岭",
  "CRO_水稻_句容",
  "GRA_高寒草甸_若尔盖",
  "GRA_人工垂穗披碱草_三江源",
]

if env_bool("BEPS_AUTORUN", true)
  Process(SITES)
  HAS_SITE_FILTER || Process(SITES_poor)
end
