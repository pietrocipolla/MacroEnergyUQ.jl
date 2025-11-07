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
@info "Configuring RCall..."
set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr)

@info "R setup completed successfully!"

@info "Installing required R packages..."

# Setup the R environment for RCall
CondaPkg.activate!(ENV)
Pkg.build("RCall")

# Install necessary R packages
R"""
    if (!require("gsaot")) {
        install.packages("gsaot", repos="https://cloud.r-project.org")
    }
"""

@info "R packages installed successfully!"