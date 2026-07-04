# NEWS

CMFD 逐小时驱动在 ChinaFlux 通量站的迭代记录（2026-06 — 2026-07）。

## NSE_CMFD_1h 系列版本对照

`case02_ChinaFlux_CMFD.jl` 一共记录了 **4 个版本**（V1/V2/V3/V4），对应输出目录：

```
OUTPUT/ALL/Bonan/NSE_CMFD_1h      # V1
OUTPUT/ALL/Bonan/NSE_CMFD_1h_V2   # V2
OUTPUT/ALL/Bonan/NSE_CMFD_1h_V3   # V3
OUTPUT/ALL/Bonan/NSE_CMFD_1h_V4   # V4
```

### 核心差别

| 版本 | 输出目录 | 站点范围 | `z_wind`（风速参考高度） | 风速来源 | `Rln_in`（长波入射） | 触发原因 / 关键改动 |
|------|---------|---------|--------------------------|---------|---------------------|---------------------|
| V1 | `NSE_CMFD_1h` | 全部 34 站 | `z_Uz`（按元数据，塔高） | 直接用 CMFD 10 m 风（视作塔高量级） | CMFD 原始值 | 基线版：与原站点驱动做对照；发现 9 站森林 ET 严重退化 |
| V2 | `NSE_CMFD_1h_V2` | 全部 34 站 | `z_Uz` | 对数廓线把 10 m 风外推到塔高 | CMFD 原始值 | §6.6 实验：验证"风速高度错配致森林退化"——结果被证伪；高冠层（h > 11.4 m）9 站直接报 `DomainError` |
| V3 | `NSE_CMFD_1h_V3` | 仅 h > 10 m 的森林站 | `z_Uz` + `safe_wind_ref` 自动抬升 | CMFD 10 m 风 | CMFD 原始值 | 用 `safe_wind_ref` 解决 V2 的 `DomainError`，对高大冠层单源对数廓线做数学兼容性修补 |
| V4 | `NSE_CMFD_1h_V4` | 仅 h < 10 m 的矮冠层站 | `10.0` m | 直接用 CMFD 10 m 风 | 强制 NaN（模型用 `cal_Rln(ϵ_air, Tair)` 估算） | §6.6 结论：矮植被应按 CMFD 实际高度设 10 m；同时验证"用气温估算长波辐射"的可行性 |

### 关键结论

- **V1→V2**：实验 §6.6 反直觉。订正风速高度并不能修复森林 ET 退化，反而变差；9 站高大森林因 `d = 0.8 h > 10 m` 导致对数廓线自变量为负，运行直接报错。
- **V2→V3**：`safe_wind_ref` 自动抬升参考高度，规避 `DomainError`，仅用于高冠层站点。
- **V2→V4**：矮植被站统一按 CMFD 实际高度设 `z_wind = 10.0` m；并放弃 CMFD 的 `Rln_in`，改由 `cal_Rln(ϵ_air, Tair)` 从气温估算长波入射。
- 全网层面 V2 与 V1 基本持平：Flux_NSE 中位数 0.525 → 0.516，等级分布由 3 优 / 16 良 / 10 一般 / 5 差 微调为 3 优 / 16 良 / 9 一般 / 6 差。
- 当前 `case02_ChinaFlux_CMFD.jl` 默认 `outdir = ".../NSE_CMFD_1h_V4"`；不存在 V5 及以后。

### 代码引用

- `case02_ChinaFlux_CMFD.jl:17` — 当前默认输出目录为 V4
- `compare_V3.jl:3` — 对照 V1/V2/V3
- `compare_V4.jl:21-22` — 对照 V2/V4
- `Plan/Report_China_FluxALL.typ` §6.5–§6.7 — 完整迭代推导与判断