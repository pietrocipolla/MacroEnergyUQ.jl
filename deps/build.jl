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
# Force override if the user explicitly requests it OR if paths don't match current preferences
force_override = get(ENV, "MACROENERGYUQ_FORCE_RHOME", "0") == "1"

if isdir(target_rhome) && isfile(target_libr)
    # Always try to set preferences, force if needed
    current_rhome = try
        load_preference(RCALL_UUID, "Rhome", nothing)
    catch
        nothing
    end
    
    # Force override if current preference points to a different R
    if current_rhome !== nothing && current_rhome != target_rhome
        @info "Detected different R installation. Forcing CondaPkg R..."
        force_override = true
    end
    
    try
        set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr; force=force_override)
        @info "RCall preferences set to Rhome=$(target_rhome) and libR=$(target_libr)."
    catch e
        if e isa ArgumentError
            @warn "RCall preferences already set to a different value. Set MACROENERGYUQ_FORCE_RHOME=1 to override, or run: julia -e 'using Preferences, UUIDs; delete_preferences!(UUID(\"6f49c342-dc21-5d91-9882-a32aef131414\"), \"Rhome\"; force=true); delete_preferences!(UUID(\"6f49c342-dc21-5d91-9882-a32aef131414\"), \"libR\"; force=true)'"
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

# Force RCall to rebuild with the CondaPkg R
@info "Building RCall with CondaPkg R environment..."
try
    Pkg.build("RCall")
catch e
    @warn "RCall build failed, trying again after clearing preferences..."
    # Delete RCall preferences and try again
    delete_preferences!(RCALL_UUID, "Rhome"; force=true)
    delete_preferences!(RCALL_UUID, "libR"; force=true)
    # Reset preferences with CondaPkg paths
    if isdir(target_rhome) && isfile(target_libr)
        set_preferences!(RCALL_UUID, "Rhome" => target_rhome, "libR" => target_libr; force=true)
    end
    Pkg.build("RCall")
end

# Install necessary R packages
R"""
    if (!require("gsaot")) {
        install.packages("gsaot", repos="https://cloud.r-project.org")
    }
"""

@info "R packages installed successfully!"