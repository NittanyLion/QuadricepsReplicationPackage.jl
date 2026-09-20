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
    pair_floor(d, p) -> Int

The counting floor of the ±pair ansatz: the smallest node count a centrally symmetric rule of
odd degree `p` in `d` dimensions can have if its even-degree moment conditions are independent.
A ±pair matches every moment of odd total degree identically and carries `d+1` unknowns for its
2 nodes, so with `M` the number of monomials of even total degree at most `p`,

    N ≥ min(2⌈M/(d+1)⌉, 2⌈(M−1)/(d+1)⌉ + 1)

the second term being the variant with one node at the center.  This is a parameter count for
that ansatz, not a lower bound for the problem: a rule whose extra symmetry makes those moment
conditions dependent can sit below it, as the Gaussian rules of the paper's first table do.  It
is printed as the `pair` column of the uniform-weight table.
"""
function pair_floor(d::Int, p::Int)
    isodd(p) || throw(ArgumentError("the pair floor as stated is for odd degree"))
    M = sum(binomial(k + d - 1, d - 1) for k in 0:2:p)
    min(2 * cld(M, d + 1), 2 * cld(M - 1, d + 1) + 1)
end

"""
    rho(n, d, p)

The quality measure of the tables: ``ρ = N^{1/d} / q`` with ``q = (p+1)/2``.  The Gauss product
grid has ``ρ = 1``; smaller is better.
"""
rho(n::Int, d::Int, p::Int) = n^(1 / d) / ((p + 1) ÷ 2)
