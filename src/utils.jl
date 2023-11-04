# function used to simplify target allocation (by dividing by GCD)
simplify(target::Vector{<:Number}) = @match target begin
    w::Vector{Int64} || w::Vector{Rational{Int64}} => Int.(w ./ gcd(w))
    _                                              => target
end


# function represents a `Number` as a string
as_string(param::Number) = @match param begin
    p::Rational{Int64} => "$(numerator(p))/$(denominator(p))"
    p::Float64         => "$(round(p, digits = 2))" 
    p::Int64           => "$(p)"
end


"""Function sets a label for a randomization procedure.

An example:
```julia-repl
julia> w = [1, 2, 3, 4]
4-element Vector{Int64}:
 1
 2
 3
 4

julia> dlr_proc = DLR(w, 2)
DLR(2): restricted randomization procedure, targeting 1:2:3:4 allocation in 4-arm trial.

julia> label = set_label(dlr_proc)
"DLR(2)"
```
"""
set_label(rnd::Union{CompleteRandomization, RestrictedRandomization}) = @match rnd begin
    r::CRD    => "CRD"
    r::TBD    => "TBD"
    r::PBD    => "PBD($(as_string(r.param)))"
    r::RAND   => "Rand"
    r::EBCD   => "EBCD($(as_string(r.p)))"
    r::ABCD   => "ABCD($(as_string(r.a)))"
    #r::GBCD   => "GBCD($(as_string(r.γ)))"
    #r::BSD    => "BSD($(as_string(r.mti)))"
    #r::BCDWIT => "BCDWIT($(as_string(r.p)), $(as_string(r.mti)))"
    #r::BUD    => "BUD($(as_string(r.λ)))"
    #r::MWUD   => "MWUD($(as_string(r.α)))"
    #r::DLR    => "DLR($(as_string(r.a)))"
    #r::DBCD   => "DBCD($(as_string(r.γ)))"
    #r::MinQD  => "MinQD($(as_string(r.η)))"
    #r::MaxEnt => "MaxEnt($(as_string(r.η)))"
end


# function used to display randomization procedure.
function Base.show(io::IO, rnd::Union{CompleteRandomization, RestrictedRandomization})
    # it is assumed that any object of `CRD` or `RestrictedRandomization` subtype
    # has the following field: `target`
    target = rnd.target
    label = set_label(rnd)

    tar = join(target, ":")  # String representation of the target allocation ratio
    ntrt = length(target)    # Calculating number of treatment arms
    class = @match rnd begin # Class of randomization: complete or restricted
        r::CRD => "complete"
        _      => "restricted"
    end

    description = "$(label): $(class) randomization procedure, targeting $(tar) allocation in $(ntrt)-arm trial."
    println(io, description)
end

