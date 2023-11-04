"""A type of rectricted randomization, representing _**A**_djustable _**B**_iased _**C**_oin _**D**_esign (_**ABCD**_).
"""
struct ABCD <: RestrictedRandomization
    target::Vector{Int64}
    a::Number

    function ABCD(target::Vector{Int64}, a::Number)
        # getting number of treatments
        ntrt = length(target)
        
        @assert ntrt == 2 "The procedure isn't implemented for multi-arm trials";
        @assert allequal(target) "The procedure isn't implemented for unequal allocation";
        
        return new(target, a)
    end
end
ABCD(a::Number) = ABCD([1, 1], a)