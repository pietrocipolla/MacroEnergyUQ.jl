@info "Setting up R environment for MacroEnergyUQ..."

using Libdl
using CondaPkg
using Preferences
using UUIDs

const RCALL_UUID = UUID("6f49c342-dc21-5d91-9882-a32aef131414")

CondaPkg.add("r")
target_rhome = joinpath(CondaPkg.envdir(), "lib", "R")
if Sys.iswindows()
    target_libr = joinpath(target_rhome, "bin", Sys.WORD_SIZE==64 ? "x64" : "i386", "R.dll")
else
    target_libr = joinpath(target_rhome, "lib", "libR.$(Libdl.dlext)")
end

try
    set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr; force=true)
    @info "RCall preferences set successfully"
catch e
    @warn "Could not set RCall preferences: $e"
end

CondaPkg.activate!(ENV)

# Now it's safe to load RCall
@info "Loading RCall..."
using RCall

# Install necessary R packages
@info "Installing required R packages..."
R"""
    if (!require("gsaot")) {
        install.packages("gsaot", repos="http://cran.rstudio.com/")
    }
"""

@info "R packages installed successfully!"

@info "✓ MacroEnergyUQ setup completed successfully!"