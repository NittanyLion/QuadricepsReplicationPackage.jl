# BigFloat verification of a quadrature rule: the worst relative (backward) moment error over
# all monomials of total degree ≤ p, and weight diagnostics.  This is the kernel that produced
# the "rel. err." column of the paper's tables and the `err_float64` / `err_extended` columns of
# the deposits, unchanged.
#
# In-place MPFR (mpfr_mul / mpfr_fma / mpfr_set) into buffers we own: the same arithmetic as
# `a[i] *= b[i]` but without allocating a fresh BigFloat per multiply.  Every buffer is built at
# an explicit precision, so the result does not depend on the task-local default precision of
# whichever thread runs the branch, and the maximum over branches does not depend on the
# number of threads.
const RNDN = Int32(0)                                     # MPFR_RNDN
mul_into!(z::BigFloat, x::BigFloat, y::BigFloat) =
    (ccall((:mpfr_mul, :libmpfr), Int32,
           (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, y, RNDN); z)
fma_into!(z::BigFloat, x::BigFloat, y::BigFloat) =                  # z += x*y
    (ccall((:mpfr_fma, :libmpfr), Int32,
           (Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Ref{BigFloat}, Int32),
           z, x, y, z, RNDN); z)
set_into!(z::BigFloat, x::BigFloat) =
    (ccall((:mpfr_set, :libmpfr), Int32,
           (Ref{BigFloat}, Ref{BigFloat}, Int32), z, x, RNDN); z)
newbuf(n, prec) = [BigFloat(0; precision = prec) for _ in 1:n]

"""
    moment1(family, e, prec)

Exact one-dimensional moment ``∫ x^e ω(x) dx`` of the family's normalized weight: the standard
normal density for `"hermite"` (``(e-1)!!`` for even `e`, 0 for odd `e`), the uniform density on
``[0,1]`` for `"legendre"` (``1/(e+1)``).
"""
function moment1(fam::AbstractString, e::Int, prec::Int)
    if fam == "hermite"
        isodd(e) && return BigFloat(0; precision = prec)
        e == 0 && return BigFloat(1; precision = prec)
        return BigFloat(prod(BigInt.(1:2:(e - 1))); precision = prec)   # (e-1)!!
    elseif fam == "legendre"
        return BigFloat(1; precision = prec) / BigFloat(e + 1; precision = prec)
    end
    throw(ArgumentError("family must be \"hermite\" or \"legendre\", not \"$fam\""))
end

# Depth-first walk over every multi-index of total degree ≤ p.  Level k holds the running
# prefix product Π_{j<k} x_j^{a_j} for every node; stepping a_k up by one costs ONE multiply per
# node, so the whole sweep costs about (#monomials × n) multiplies instead of (#monomials × n × d).
function walk!(k::Int, rem::Int, src::Vector{BigFloat}, srcabs::Vector{BigFloat},
               W::Vector{Vector{BigFloat}}, WA::Vector{Vector{BigFloat}},
               expo::Vector{Int}, Xc::Vector{Vector{BigFloat}},
               AXc::Vector{Vector{BigFloat}}, w::Vector{BigFloat},
               aw::Vector{BigFloat}, m1::Vector{BigFloat}, d::Int,
               acc::Vector{BigFloat}, prec::Int)
    if k > d
        q, s, ex, tmp, zer, one_ = acc[1], acc[2], acc[3], acc[4], acc[5], acc[6]
        set_into!(q, zer)
        set_into!(s, zer)
        @inbounds for t in eachindex(w)
            fma_into!(q, w[t], src[t])
            fma_into!(s, aw[t], srcabs[t])
        end
        set_into!(ex, m1[expo[1] + 1])
        @inbounds for j in 2:d
            mul_into!(tmp, ex, m1[expo[j] + 1]); set_into!(ex, tmp)
        end
        return Float64(abs(q - ex) / max(s, one_))
    end
    worst = 0.0
    Wk, WAk = W[k], WA[k]
    @inbounds for t in eachindex(Wk)
        set_into!(Wk[t], src[t]); set_into!(WAk[t], srcabs[t])
    end
    for e in 0:rem
        expo[k] = e
        worst = max(worst, walk!(k + 1, rem - e, Wk, WAk, W, WA, expo, Xc, AXc,
                                 w, aw, m1, d, acc, prec))
        if e < rem
            xk, axk = Xc[k], AXc[k]
            @inbounds for t in eachindex(Wk)
                mul_into!(Wk[t], Wk[t], xk[t]); mul_into!(WAk[t], WAk[t], axk[t])
            end
        end
    end
    expo[k] = 0
    worst
end

"""
    verify_rule(nodes, weights, p, family; precision = 192) -> NamedTuple

Verify a rule in BigFloat arithmetic at `precision` bits.  `nodes` is an `n × d` matrix and
`weights` a vector of length `n`, both `BigFloat`; `family` is `"hermite"` (standard normal
weight on ``ℝ^d``) or `"legendre"` (uniform weight on ``[0,1]^d``).

The result has the fields

* `err`: ``\\max_{|a| ≤ p} |∑_s w_s x_s^a − E\\,x^a| / \\max(∑_s |w_s|\\,|x_s^a|, 1)``, the largest
  relative error over **all** monomials of total degree at most `p`;
* `minw`, `maxw`: the smallest and largest weight;
* `sumw`, `sumwdev`: the sum of the weights and its distance from 1;
* `minnode`, `maxnode`, `maxabs`: the range of the node coordinates.
"""
function verify_rule(X::Matrix{BigFloat}, w::Vector{BigFloat}, p::Int,
                     fam::AbstractString; precision::Int = 192)
    prec = precision
    n, d = size(X)
    n == length(w) || throw(DimensionMismatch("$n nodes but $(length(w)) weights"))
    Xc  = [[X[t, k] for t in 1:n] for k in 1:d]
    AXc = [[abs(X[t, k]) for t in 1:n] for k in 1:d]
    aw  = [abs(wt) for wt in w]
    m1  = [moment1(fam, e, prec) for e in 0:p]
    parts = zeros(Float64, p + 1)
    Threads.@threads for e1 in 0:p
        W  = [newbuf(n, prec) for _ in 1:d]
        WA = [newbuf(n, prec) for _ in 1:d]
        acc = newbuf(6, prec)          # q, s, exact, tmp, 0, 1
        expo = zeros(Int, d)
        src  = newbuf(n, prec); srcabs = newbuf(n, prec)
        one_ = BigFloat(1; precision = prec)
        set_into!(acc[6], one_)         # acc[5] stays zero
        for t in 1:n
            set_into!(src[t], one_); set_into!(srcabs[t], one_)
        end
        # apply x_1^e1 to the level-1 prefix, then walk dimensions 2..d
        for _ in 1:e1, t in 1:n
            mul_into!(src[t], src[t], Xc[1][t]); mul_into!(srcabs[t], srcabs[t], AXc[1][t])
        end
        expo[1] = e1
        parts[e1 + 1] = d == 1 ?
            walk!(2, 0, src, srcabs, W, WA, expo, Xc, AXc, w, aw, m1, d, acc, prec) :
            walk!(2, p - e1, src, srcabs, W, WA, expo, Xc, AXc, w, aw, m1, d, acc, prec)
    end
    sumw = sum(w)
    (err     = maximum(parts),
     minw    = Float64(minimum(w)),
     maxw    = Float64(maximum(w)),
     sumw    = Float64(sumw),
     sumwdev = Float64(abs(sumw - 1)),
     minnode = Float64(minimum(X)),
     maxnode = Float64(maximum(X)),
     maxabs  = Float64(maximum(abs, X)))
end
