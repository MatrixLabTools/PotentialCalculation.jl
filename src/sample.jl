module sample

export line_sampler, adaptive_line_sampler, sample_multiple_adaptive_lines

using LinearAlgebra
using ..calculators, ..clusters, ..atoms, ..unitconversions

# Energy conversion. TODO make this use uniconv from PotentialFit
#function energy_conv(e, unit)
#    if unit in ["cm-1", "cm^-1"]
#        return e/(27.2113834*8065.54477)
#    elseif unit in ["ev", "eV"]
#        return e/27.2113834
#    else
#        return e
#    end
#end


function random_rotation!(c::AbstractCluster)
    theta = 2*pi*rand(3)
    rotate_x!(c,theta[1])
    rotate_y!(c,theta[2])
    rotate_z!(c,theta[3])
end

function random_rotation_matrix()
    theta = 2*π*rand(3)
    Rx = Matrix{Float64}(I, 3,3)
    Rx[2,2] = cos(theta[1])
    Rx[3,3] = Rx[2,2]
    Rx[3,2] = -sin(theta[1])
    Rx[2,3] = -Rx[3,2]

    Ry = Matrix{Float64}(I, 3,3)
    Ry[1,1] = cos(theta[2])
    Ry[3,3] = Ry[1,1]
    Ry[3,1] = sin(theta[2])
    Ry[1,3] = -Ry[3,1]

    Rz = Matrix{Float64}(I, 3,3)
    Rz[1,1] = cos(theta[3])
    Rz[2,2] = Rz[1,1]
    Rz[2,1] = -sin(theta[3])
    Rz[1,2] = -Rz[2,1]

    return Rx * Ry * Rz
end


function _line_sampler(cl1::Cluster, cl2::Cluster, u, step, npoints)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)
    out1 = Array{typeof(c1)}(undef, npoints)
    out2 = Array{typeof(c2)}(undef, npoints)
    for i in eachindex(out1)
        out1[i] = deepcopy(c1)
        out2[i] = deepcopy(c2)
        move!(c2, step*u)
    end
    return out1, out2
end

function line_sampler(cl1::Cluster, cl2::Cluster; npoints=10, mindis=2.0, maxdis=9.0)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)

    center_cluster!(c1)
    center_cluster!(c2)

    step = (maxdis - mindis) / npoints
    u = random_rotation_matrix() * [1.0, 0.0, 0.0]
    move!(c2, mindis * u )
    while minimum(distances(c1,c2)) < mindis
        move!(c2, step * u)
    end
     return _line_sampler(c1,c2,u,step,npoints)
end


"""
  adaptive_line_sampler(cal::Calculator, cl1::Cluster, cl2::Cluster, max_e=0; unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0)

Calculates potential on a line type distances that is suitable for visualization.
Closest distance used in calculations is calculated adaptively.

# Arguments
- `max_e` : Point that is closest and has less energy than this will be starting point for a line
- `unit` : Unit in which `max_e` is given
- `npoints` : Number of points in potential
- `maxdis` : Maximum distance in potential calculation
- `sstep` : Search step size for adaptive search
- `startdistance` : Distance from which adaptive seach is started
"""
function adaptive_line_sampler(cal::Calculator, cl1::Cluster, cl2::Cluster, max_e=0;
                               unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0)
    @info "Starting adaptive_line_sampler"
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)

    emax = energy_from(max_e, unit)
    center_cluster!(c1)
    center_cluster!(c2)

    u = random_rotation_matrix() * [1.0, 0.0, 0.0]
    r = 0.5*startdistance
    move!(c2, 0.5*startdistance*u)
    while minimum(distances(c1,c2)) < startdistance
        move!(c2, sstep * u)
        r += sstep
    end

    e = bsse_corrected_energy(cal, [c1], [c2]) .- emax
    ctemp = [deepcopy(c2)]
    fless = false
    fmore = false
    if e[1] > 0
        fmore = true
    else
        fless = true
    end
    while fmore == false || fless == false
        @debug "Step energy is $(e[end])"
        @debug "Mindistance is $(minimum(distances(c1,c2)))"
        if e[end] < 0
            move!(c2, -1*sstep*u)
            r -= sstep
        else
            move!(c2, sstep*u)
            r += sstep
        end
        tmp = bsse_corrected_energy(cal, [c1], [c2]) .- emax
        push!(e,  tmp...)
        push!(ctemp, deepcopy(c2))
        if e[end] > 0
            fmore = true
        else
            fless = true
        end
    end
    if length(e) == 1
        error("Search fail!")
    end
    @info "adaptive_line_sampler found initial point in $(length(e)) steps"
    @debug "e $(e)"
    out2 = [ctemp[e.<0][end]]
    eout = [e[e.<0][end] + emax]
    out1 = [c1]
    c2 = deepcopy(out2[1])
    @debug "Distances  $(distances(c1,c2))"
    step = (maxdis - r ) / npoints
    move!(c2,step*u)
    tmp = _line_sampler(c1, c2[1], u, step, npoints-1)
    et = bsse_corrected_energy(cal, tmp[1], tmp[2])
    push!(eout, et...)
    push!(out1, tmp[1]...)
    push!(out2, tmp[2]...)
    @debug minimum.(distances.(out1,out2))
    return Dict("Energy" => eout, "Points" =>  out1 .+ out2,  "Mindis" =>  minimum(distances(out1[1],out2[1])))
end


function sample_multiple_adaptive_lines(cal::Calculator, cl1::Cluster, cl2::Cluster, nlines, max_e=0;
                               unit="cm-1", npoints=10, maxdis=9.0, sstep=0.1, startdistance=3.0)
    c1 = deepcopy(cl1)
    c2 = deepcopy(cl2)
    out = Vector{Dict}(undef, nlines)
    rtmp = Float64[]
    sr = startdistance
    for i in 1:nlines
        tmp = adaptive_line_sampler(cal, c1, c2, max_e; unit=unit, npoints=10, maxdis=9.0, startdistance=sr)
        push!(rtmp, tmp["Mindis"])
        sr = sum(rtmp) / length(rtmp)
        out[i] = tmp
    end
    return out
end


end  # module sample
