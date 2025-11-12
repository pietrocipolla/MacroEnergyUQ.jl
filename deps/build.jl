using Libdl
using CondaPkg
using Preferences
using UUIDs
import Pkg

const RCALL_UUID = UUID("6f49c342-dc21-5d91-9882-a32aef131414")

@info "Setting up R environment for MacroEnergyUQ..."

# Install R through conda
CondaPkg.add("r")

# Activate conda environment
CondaPkg.activate!(ENV)

# Set up R paths from CondaPkg
target_rhome = joinpath(CondaPkg.envdir(), "lib", "R")
if Sys.iswindows()
    target_libr = joinpath(target_rhome, "bin", Sys.WORD_SIZE==64 ? "x64" : "i386", "R.dll")
else
    target_libr = joinpath(target_rhome, "lib", "libR.$(Libdl.dlext)")
end

# Clear any existing RCall preferences that might point to system R
@info "Clearing existing RCall preferences..."
try
    delete_preferences!(RCALL_UUID, "Rhome"; force=true)
    delete_preferences!(RCALL_UUID, "libR"; force=true)
    @info "Cleared existing RCall preferences"
catch e
    @info "No existing preferences to clear"
end

# Set RCall preferences to use CondaPkg R
@info "Configuring RCall to use CondaPkg R..."
if isdir(target_rhome) && isfile(target_libr)
    set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr; force=true)
    @info "RCall preferences set to:"
    @info "  Rhome = $(target_rhome)"
    @info "  libR  = $(target_libr)"
else
    error("CondaPkg R installation not found at $(target_rhome). Please ensure CondaPkg.add(\"r\") succeeded.")
end

# Build RCall with the correct preferences
@info "Building RCall..."
Pkg.build("RCall")

# Now it's safe to load RCall
@info "Loading RCall..."
using RCall

# Verify RCall is using the correct R
rcall_rhome = rcopy(R"R.home()")
if rcall_rhome != target_rhome
    @warn "RCall is using R from $(rcall_rhome) instead of $(target_rhome)"
    @warn "You may need to restart Julia and rebuild: julia -e 'using Pkg; Pkg.build(\"MacroEnergyUQ\")'"
else
    @info "✓ RCall successfully configured to use CondaPkg R"
end

# Install necessary R packages
@info "Installing required R packages..."
R"""
    if (!require("gsaot")) {
        install.packages("gsaot", repos="https://cloud.r-project.org")
    }
"""

@info "R packages installed successfully!"

@info "✓ MacroEnergyUQ setup completed successfully!"