"""
    moller_bound(d, p) -> Int

Möller's lower bound on the number of nodes of a cubature rule of odd degree `p` in `d`
dimensions for a centrally symmetric weight (both weights of the paper qualify).  It binds every
rule, symmetric or not, so a rule that attains it is proven minimal.

H. M. Möller, *Lower bounds for the number of nodes in cubature formulae*, in Numerische
Integration, ISNM 45, Birkhäuser 1979, in the form stated by Orive, Santos-León and Spalević,
ETNA 53 (2020), Theorem 1.1: with ``p = 2s − 1``,

    N ≥ C(d+s−1, d) + Σ_{k=1}^{d−1} 2^{k−d} C(k+s−1, k)            (s even)
    N ≥ C(d+s−1, d) + Σ_{k=1}^{d−1} (1 − 2^{k−d}) C(k+s−2, k)      (s odd)

evaluated in exact rational arithmetic and rounded up.
"""
function moller_bound(d::Int, p::Int)
    isodd(p) || throw(ArgumentError("Möller's bound as stated is for odd degree"))
    s = (p + 1) ÷ 2
    tot = Rational{BigInt}(binomial(big(d + s - 1), big(d)))
    half = Rational{BigInt}(1, 2)
    for k in 1:(d - 1)
        if iseven(s)
            tot += half^(d - k) * binomial(big(k + s - 1), big(k))
        else
            tot += (1 - half^(d - k)) * binomial(big(k + s - 2), big(k))
        end
    end
    Int(ceil(tot))
end

"""
    rho(n, d, p)

The quality measure of the tables: ``ρ = N^{1/d} / q`` with ``q = (p+1)/2``.  The Gauss product
grid has ``ρ = 1``; smaller is better.
"""
rho(n::Int, d::Int, p::Int) = n^(1 / d) / ((p + 1) ÷ 2)
