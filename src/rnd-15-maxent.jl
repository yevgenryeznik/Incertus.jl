"""A type of rectricted randomization, representing _**Max**_imum _**Ent**_opy Constraint Balance Randomization (_**MaxEnt**_).

`MaxEnt(η)` command initializes _MaxEnt design_ with a parameter equal to `η`,
targeting `1:1` allocation.

`MaxEnt(w, η)` command initializes _MaxEnt design_ with a parameter equal to `η`,
targeting allocation specified by `w`.

An output of the command is an isntance of MaxEnt.
"""
struct MaxEnt <: RestrictedRandomization 
    target::Vector{Int64}
    parameter::Number

    function MaxEnt(target::Vector{Int64}, parameter::Number)
        @assert 0 <= parameter <= 1 "The procedure's parameter, `η`, must belong to [0; 1]";
  
        return new(target, parameter)
    end
end
MaxEnt(parameter::Number) = MaxEnt([1, 1], parameter)


"""Function calculates allocation probabilities for MaxEnt, given treatment numbers.
# Call
- `allocation_prb(rnd, N)`

# Arguments
- `rnd::MaxEnt`: an object, representing Maximum Entropy Constraint Balanced Randomiztion.
- `N::Vector{Int64}`: a vector of current treatment numbers.
"""
function allocation_prb(rnd::MaxEnt, N::Vector{Int64})
    # target allocation
    w = rnd.target
    
    # parameter of the randomization procedure (MWUD)
    η = rnd.parameter

    # getting number of treatment arms
    ntrt = length(w)

    # current subject's ID
    j = sum(N) + 1

    # target allocation proportions
    ρ = w ./ sum(w)

    # the hypothetical "lack of balance"
    B = zeros(ntrt)
    N1 = N
    for k ∈ 1:ntrt
        N1[k] += 1
        B[k] = maximum(abs.(N1./j - ρ))
        N1[k] -= 1 
    end

    # probabilities of tretament assignments
    if var(B) <= 1e-16
        return ρ
    else
        # settin a nonlinear function to find its zero
        function f(μ, p) 
            B, ρ, η = p
            return minimum(B)*η + (1-η)*sum(ρ.*B) - sum(B.*ρ.*exp.(-μ.*B))/sum(ρ.*exp.(-μ.*B))
        end
        # setting and solving a zero-fincding problem to find μ
        zfp = ZeroProblem(x -> f(x, (B, ρ, η)), 0)
        μ = solve(zfp)

        p_num = ρ .* exp.(-μ .* B)
        
        return p_num/sum(p_num)
    end
end