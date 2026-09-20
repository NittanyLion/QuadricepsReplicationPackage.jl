"""
    QuadricepsReplicationPackage

Replication package for J. Pinkse, *Positive weight Hermite and Legendre quadrature rules*
(2026).  The paper's contribution is two tables of positive-weight quadrature rules, for the
Gaussian weight ``N(0, I_d)`` and for the uniform weight on ``[0,1]^d``, in ``d = 2, …, 5``
dimensions.  This package holds every rule of those tables in double precision and rebuilds
the tables from them:

```julia
using QuadricepsReplicationPackage
results = replicate()        # all cells; minutes, threaded
replicated(results)          # true when every check passed
```

See [`replicate`](@ref), [`check_cell`](@ref), [`check_extended`](@ref), [`check_float128`](@ref), [`verify_rule`](@ref)
and [`moller_bound`](@ref).
"""
module QuadricepsReplicationPackage

using Printf, SHA

export replicate, replicated, check_cell, check_extended, check_float128, verify_rule, moller_bound, pair_floor, rho,
       cells, load_rule, load_rule64, paper_rows, table_row, datadir, Cell, CellResult, PaperRow

include("verify.jl")
include("moller.jl")
include("rules.jl")
include("tables.jl")
include("replicate.jl")

end # module
