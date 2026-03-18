# Fallback definitions used when the optional RCall extension is not loaded.
# Real implementations are provided by ext/MacroEnergyUQRCallExt.jl.

const RCALL_MISSING_MSG = "RCall is not available. Install and load RCall to enable R-based features.\n" *
                         "Example: ] add RCall"

function ot_indices(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.ot_indices(args...; kwargs...)
    end
    throw(ArgumentError(RCALL_MISSING_MSG))
end

function ot_indices_wb(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.ot_indices_wb(args...; kwargs...)
    end
    throw(ArgumentError(RCALL_MISSING_MSG))
end

function ot_indices_1d(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.ot_indices_1d(args...; kwargs...)
    end
    throw(ArgumentError(RCALL_MISSING_MSG))
end

function transport_ot(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.transport_ot(args...; kwargs...)
    end

    throw(ArgumentError(RCALL_MISSING_MSG))
end

function transport_ranking(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.transport_ranking(args...; kwargs...)
    end

    throw(ArgumentError(RCALL_MISSING_MSG))
end

function irrelevance_threshold(args...; kwargs...)
    ext = Base.get_extension(@__MODULE__, :MacroEnergyUQRCallExt)
    if !isnothing(ext)
        return ext.irrelevance_threshold(args...; kwargs...)
    end

    throw(ArgumentError(RCALL_MISSING_MSG))
end
