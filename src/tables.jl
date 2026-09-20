"""
    PaperRow

One row of a table in the paper, as typeset: `N` nodes, `rho`, `err` (the "rel. err." column as
printed), `moller` (Möller's bound), `pair` (the ±pair counting floor, printed in the uniform-weight
table only and `nothing` in the Gaussian one), `prev` (the best count in the literature,
`nothing` where only the product grid exists), `src` (the code of the source of `prev`, starred
when the paper's rule was warm-started from theirs), and `shade` (`""`, `"light"` or `"dark"`).
"""
struct PaperRow
    family::String
    d::Int
    p::Int
    N::Int
    rho::String
    err::String
    moller::Union{Int,Nothing}
    pair::Union{Int,Nothing}
    prev::Union{Int,Nothing}
    src::String
    shade::String
end

const SHADE = Dict("" => "", "light" => "\\rowcolor{gray!15}", "dark" => "\\rowcolor{gray!40}")

"""
    paper_rows(family; dir = datadir()) -> Dict{Tuple{Int,Int},PaperRow}

The rows of the paper's table for `family` (`"gh"` or `"le"`), keyed by `(d, p)`, read from the
LaTeX source of the table as it appears in the paper (`data/paper/best_known_<family>.tex`).
"""
function paper_rows(fam::AbstractString; dir::AbstractString = datadir())
    out = Dict{Tuple{Int,Int},PaperRow}()
    d = 0
    # the uniform-weight table carries one column more than the Gaussian one: the pair floor,
    # between Möller's bound and prev (2026-09-19)
    for l in eachline(joinpath(dir, "paper", "best_known_$fam.tex"))
        md = match(r"^\\multicolumn\{[78]\}\{@\{\}l\}\{\$d = (\d+)\$", l)
        md === nothing || (d = parse(Int, md[1]); continue)
        m = match(r"^(\\rowcolor\{gray!(\d+)\})?\s*(\d+) & (\d+) & ([\d.]+) & (\S+) & (\S+) & (\S+)(?: & (\S+))? & (.*?) \\\\\s*$", l)
        m === nothing && continue
        shade = m[2] === nothing ? "" : m[2] == "40" ? "dark" : "light"
        num(s) = s === nothing || s == "---" ? nothing : parse(Int, s)
        p = parse(Int, m[3])
        pair, prev = m[9] === nothing ? (nothing, num(m[8])) : (num(m[8]), num(m[9]))
        out[(d, p)] = PaperRow(fam, d, p, parse(Int, m[4]), m[5], m[6], num(m[7]), pair, prev,
                               String(m[10]), shade)
    end
    out
end

texerr(e) = e == 0 ? "0" : @sprintf("%.1e", e)

"""
    table_row(cell, err, prev, src) -> String

The LaTeX row of the paper's table for `cell`, rebuilt from the rule: `N` is the number of data
lines of the rule file, `ρ` and Möller's bound are computed, `err` is the verified error, and
the shading follows from `N`, the bound and `prev`.  `prev` and `src` describe the literature
and are taken as given.
"""
function table_row(c::Cell, err::Float64, prev::Union{Int,Nothing}, src::AbstractString)
    mb = moller_bound(c.d, c.p)
    shade = c.n == mb ? "dark" : (prev !== nothing && c.n ≥ prev) ? "light" : ""
    pair = c.family == "le" ? string(pair_floor(c.d, c.p), " & ") : ""      # uniform weight only
    string(SHADE[shade], c.p, " & ", c.n, " & ", @sprintf("%.2f", rho(c.n, c.d, c.p)), " & ", texerr(err), " & ",
           mb, " & ", pair, prev === nothing ? "---" : prev, " & ", src, " \\\\")
end

"The row as the paper typesets it, for comparison with `table_row`."
paper_row_tex(r::PaperRow) =
    string(SHADE[r.shade], r.p, " & ", r.N, " & ", r.rho, " & ", r.err, " & ",
           r.moller === nothing ? "---" : r.moller, " & ", r.pair === nothing ? "" : string(r.pair, " & "),
           r.prev === nothing ? "---" : r.prev, " & ", r.src, " \\\\")
