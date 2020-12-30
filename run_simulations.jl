# Perform actomyosin network simulations
# Alex Tam, 12/10/2020

# Load packages
using Revise # Modify code without restarting Julia
using Parameters # Tools for storing parameters in data structures
using LinearAlgebra # Matrix and vector operations
using StatsBase # Sampling without replacement
using Plots # Plotting library
using LaTeXStrings # Display LaTeX output
using ForwardDiff # Automatic differentiation tools
using Optim # Optimisation routines
using LineSearches # Line searches for optimisation
using Printf # C-style printing macros
using Distributions # Probability distributions

# Include code from external files
include("model_parameters.jl")
include("actomyosin_network.jl")
include("State.jl")
include("initial_condition.jl")
include("Actin_Filament.jl")
include("periodic.jl")
include("Cross_Link.jl")
include("Myosin_Motor.jl")
include("intersection_search.jl")
include("spatial_statistics.jl")
include("network_force.jl")
include("draw_network.jl")
include("optimise_network.jl")
include("energy.jl")

# Specify parameters
parN = Numerical_Parameters(xTol = 0, fTol = 0, nT = 201); # Initialise struct of numerical parameters
parA = Actin_Properties(nSeg = 1, LSeg = 1.0, lambda_a = 10, k = 100000000); # Initialise struct of actin filament properties
parM = Myosin_Properties(k = 100000000); # Initialise struct of myosin motor properties

# Run simulations
@time state, af, mm, xl, Force, Curvature, Dipole = actomyosin_network(parN, parA, parM);