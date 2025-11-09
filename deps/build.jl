using Libdl
using CondaPkg
using Preferences
using UUIDs
using RCall
import Pkg

const RCALL_UUID = UUID("6f49c342-dc21-5d91-9882-a32aef131414")

@info "Setting up R environment for MacroEnergyUQ..."

# Install R through conda
CondaPkg.add("r")

# Set up R paths
target_rhome = joinpath(CondaPkg.envdir(), "lib", "R")
if Sys.iswindows()
    target_libr = joinpath(target_rhome, "bin", Sys.WORD_SIZE==64 ? "x64" : "i386", "R.dll")
else
    target_libr = joinpath(target_rhome, "lib", "libR.$(Libdl.dlext)")
end

# Set RCall preferences
@info "Configuring RCall preferences..."
# Only set preferences if paths exist; don't crash if already set to a different value.
force_override = get(ENV, "MACROENERGYUQ_FORCE_RHOME", "0") == "1"
if isdir(target_rhome) && isfile(target_libr)
    try
        set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr; force=force_override)
        @info "RCall preferences set to Rhome=$(target_rhome) and libR=$(target_libr)."
    catch e
        if e isa ArgumentError
            @info "RCall preferences already set. Keeping existing values. Set MACROENERGYUQ_FORCE_RHOME=1 to override."
        else
            rethrow()
        end
    end
else
    @warn "Target R not found at $(target_rhome). Skipping setting RCall preferences."
end

@info "R setup completed successfully!"

@info "Installing required R packages..."

# Ensure the R environment is prepared for RCall build
CondaPkg.activate!(ENV)
Pkg.build("RCall")

# Install necessary R packages
R"""
    if (!require("gsaot")) {
        install.packages("gsaot", repos="https://cloud.r-project.org")
    }
"""

@info "R packages installed successfully!"