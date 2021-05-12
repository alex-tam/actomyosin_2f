# Simulate a 2D actomyosin network
# Alex Tam, 12/10/2020

# Simulate network
"Control function for actomyosin network simulations"
function actomyosin_network(parN, parA, parM)
    # Specify domain width
    Lxx::Float64 = 2.0; Lxy::Float64 = 0; Lyx::Float64 = 0; Lyy::Float64 = 2.0; # Actual domain widths
    # Pre-allocate
    Force = [[0.0, 0.0, 0.0, 0.0] for idx in 1:parN.nT]; # [pN] Network force
    Bulk_Stress = [0.0 for idx in 1:parN.nT]; # Bulk stress
    Motor_Pos = [0.0 for idx in 1:parN.nT]; # Bulk stress
    Curvature = [0.0 for idx in 1:parN.nT]; # Mean filament curvature
    Index = [0.0 for idx in 1:parN.nT]; # Two-filament index
    Theta = [0.0 for idx in 1:parN.nT]; # Angle between filaments
    # Generate initial conditions
    state = State{Float64}(Vector{Vector{Vector}}(), Vector{Vector}()); # Initialise empty State struct
    mm = Vector{Myosin_Motor}(); # Pre-allocate empty myosin motors
    af, state = actin_ic(state, parN, parA, Lxx, Lyy); # Initialise actin filaments
    xl = intersection_search(state, af, mm); # Initialise cross-links
    mm, xl, state = myosin_ic(state, mm, parN, xl, Lxx, Lxy, Lyx, Lyy); # Initialise myosin motors
    state_old = state; # Store initial state to compute force
    # Time-stepping
    animation = @animate for i = 1:parN.nT
        # Compute force and draw network
        Force[i] = network_force(state, state_old, af, xl, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy);
        draw_network(state, af, mm, parN, parA, Force[i], Lxx, Lxy, Lyx, Lyy);
        savefig("2f-$i.svg"); # Save image
        # Spatial statistics
        if parA.nSeg > 1
            Curvature[i] = curvature(af, state, Lxx, Lxy, Lyx, Lyy); # Compute curvature at current time step
        else
            Curvature[i] = 0;
        end
        Index[i], Theta[i] = two_filament_index(mm, state, Lxx, Lxy, Lyx, Lyy); # Compute two-filament index at current time step
        Bulk_Stress[i] = 0.5*(Force[i][1]/Lyy + Force[i][4]/Lxx);
        Motor_Pos[i] = 0.5*(state.mp[1][1]+state.mp[1][2]);
        # Compute solution
        if i != parN.nT # Ensure correct looping sequence
            # Turnover and polymerisation
            af = segment_translations(state, af); # Update filament translations       
            state_old = state; # Store current state for energy functional
            # Compute network solution
            new_dof = optimise_network(state_old, af, mm, parN, parA, parM, Lxx, Lxy, Lyx, Lyy);
            state = build_state(new_dof, af, mm); # Construct State from vector
        end
    end
    # Output
    gif(animation, "2f.gif", fps = 50) # Save .gif of the network
    # Plot quantities versus time
    draw_stress(parN, Force, parN.nT, Lxx, Lyy); savefig("2f_stress.svg"); # Stress components
    draw_bulk_stress(parN, Force, parN.nT, Lxx, Lyy); savefig("2f_bulk_stress.svg"); # Bulk stress
    draw_angle(parN, Theta, parN.nT); savefig("2f_angle.svg"); # Angle
    # Write data to files
    writedlm("2f_angle.csv", Theta);
    writedlm("2f_index.csv", Index);
    writedlm("2f_bulk_stress.csv", Bulk_Stress);
    writedlm("2f_motor_pos.csv", Motor_Pos);
    # Compute time-integrated statistics
    Curvature_Int, Index_Int = integrated_statistics(parN, Curvature, Index);
    return state, af, mm, Force, Curvature_Int, Index_Int
end