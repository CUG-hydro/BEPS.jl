using DocStringExtensions: TYPEDFIELDS
using DataFrames: DataFrame

include("LabelledArrays.jl")
include("helpers.jl")
include("s_coszs.jl")

include("AirLayer.jl") # TODO: test this part
# include("fill_meteo.jl")
include("lai2.jl")
include("Vcmax.jl")
include("snow_density.jl")

export s_coszs, lai2, fill_meteo!
