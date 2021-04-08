# Data structures and methods for actin filaments
# Alex Tam, 12/10/2020

"Data structure for actin filament segments"
mutable struct Segment
    index::Int # Integer segment index
    L_eq::Float64 # [μm] Segment equilibrium length
    t::Vector{Vector{Int}} # Translations required by periodic BCs
end

"Mutable data structure for actin filaments"
mutable struct Actin_Filament
    index::Int # Integer filament index
    segments::Vector{Segment} # List of segments on actin filament
    "Inner constructor to generate actin filaments and update State"
    function Actin_Filament(s, ind, parA, Lxx, Lyy)
        # 1. Generate actin filament segments
        segments = Vector{Segment}(undef, parA.nSeg); # Pre-allocate
        for i = 1:parA.nSeg
            segments[i] = Segment(i, parA.LSeg, []); # Create segments
        end
        # 2. Generate centre position and orientation
        angle = pi/4 + (ind-1)*pi/2; # Prescribe actin filament orientation
        # 3. Calculate minus end position
        nodes = Vector{}(undef, parA.nSeg+1); # Pre-allocate
        nx = 0.5; ny = 0.5; # Dimensionless actin filament minus end position
        nodes[1] = [nx, ny]; # Store minus end position
        # 4. Calculate node positions
        for i = 1:parA.nSeg
            nx += segments[i].L_eq*cos(angle)/Lxx;
            ny += segments[i].L_eq*sin(angle)/Lyy;
            nodes[i+1] = [nx, ny];
        end
        push!(s.an, nodes); # Push nodes to State
        # 5. Apply periodic BCs to each segment
        for i = 1:parA.nSeg
            segments[i].t = periodic(s.an[ind][i+1], s.an[ind][i]); # Store required translations
        end
        return new(ind, segments)
    end
end

"Obtain dimensional segment lengths"
function get_segment_lengths(f::Actin_Filament, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    Ls = Vector{T}(undef, length(f.segments)); # Pre-allocate segment lengths
    for seg in f.segments
        mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy); # Obtain segment nodes (dimensional, un-translated)
        Ls[seg.index] = sqrt((px - mx)^2 + (py - my)^2); # Compute segment length
    end
    return Ls
end

"Obtain dimensional length of a particular segment"
function get_segment_length(f::Actin_Filament, s::State{T}, seg::Segment, Lxx, Lxy, Lyx, Lyy) where {T}
    mx, my, px, py = get_segment_nodes(f, seg, s, Lxx, Lxy, Lyx, Lyy); # Obtain segment nodes (dimensional, un-translated)
    Ls = sqrt((px - mx)^2 + (py - my)^2); # Compute segment length
    return Ls
end

"Obtain dimensionless, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}) where {T}
    mx = s.an[f.index][seg.index][1];
    my = s.an[f.index][seg.index][2];
    px = s.an[f.index][seg.index+1][1];
    py = s.an[f.index][seg.index+1][2];
    return mx, my, px, py
end

"Obtain dimensionless, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}) where {T}
    mx = s.an[f.index][seg.index][1] - t[1];
    my = s.an[f.index][seg.index][2] - t[2];
    px = s.an[f.index][seg.index+1][1] - t[1];
    py = s.an[f.index][seg.index+1][2] - t[2];
    return mx, my, px, py
end

"Obtain dimensional, un-translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = s.an[f.index][seg.index][1]*Lxx + s.an[f.index][seg.index][2]*Lyx;
    my = s.an[f.index][seg.index][1]*Lxy + s.an[f.index][seg.index][2]*Lyy;
    px = s.an[f.index][seg.index+1][1]*Lxx + s.an[f.index][seg.index+1][2]*Lyx;
    py = s.an[f.index][seg.index+1][1]*Lxy + s.an[f.index][seg.index+1][2]*Lyy;
    return mx, my, px, py
end

"Obtain dimensional, translated node positions"
function get_segment_nodes(f::Actin_Filament, seg::Segment, s::State{T}, t::Vector{Int}, Lxx, Lxy, Lyx, Lyy) where {T}
    mx = (s.an[f.index][seg.index][1] - t[1])*Lxx + (s.an[f.index][seg.index][2] - t[2])*Lyx;
    my = (s.an[f.index][seg.index][1] - t[1])*Lxy + (s.an[f.index][seg.index][2] - t[2])*Lyy;
    px = (s.an[f.index][seg.index+1][1] - t[1])*Lxx + (s.an[f.index][seg.index+1][2] - t[2])*Lyx;
    py = (s.an[f.index][seg.index+1][1] - t[1])*Lxy + (s.an[f.index][seg.index+1][2] - t[2])*Lyy;
    return mx, my, px, py
end