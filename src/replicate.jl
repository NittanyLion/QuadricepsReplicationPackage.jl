const GATE = 1e-11                                   # the paper's acceptance threshold
const EXT_TARGET = Dict("gh" => 1e-68, "le" => 1e-68)  # targets of the extended-precision files (both weights ship 80 digits)
const EPS128 = 2.0^-112                              # machine epsilon of IEEE binary128 (Float128), 1.93e-34
const F128_BOUND = 10 * EPS128                       # an extended-precision file rounded to binary128 must be this exact

"""
    recorded_extended(family; dir = datadir()) -> Dict{(d, p, N) => error}

The error of each rule's extended-precision file as the deposit records it (`err_extended` in
`data/<family>/summary.csv`).  The "rel. err." column of the paper's tables is that error, and
the extended-precision files are not part of this package (they are in the Zenodo deposits), so
a table row is rebuilt from the recorded figure; `replicate(extended = …)` recomputes it from
the file and requires the two to agree.
"""
function recorded_extended(fam::AbstractString; dir::AbstractString = datadir())
    out = Dict{Tuple{Int,Int,Int},Float64}()
    f = joinpath(dir, fam, "summary.csv"); isfile(f) || return out
    L = readlines(f); h = split(L[1], ","); i = findfirst(==("err_extended"), h)
    i === nothing && return out
    back = length(h) - i                       # counted from the end: earlier quoted fields hold commas
    for l in L[2:end]
        m = match(r"^\w+,(\d+),(\d+),\d+,(\d+),", l); m === nothing && continue
        e = tryparse(Float64, split(l, ",")[end-back]); e === nothing && continue
        out[(parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3]))] = e
    end
    out
end

"""
    CellResult

What `replicate` found for one cell.  `problems` is empty when every check passed.
"""
struct CellResult
    cell::Cell
    err::Float64
    minw::Float64
    sumwdev::Float64
    row::String
    paper::Union{String,Nothing}
    err_extended::Union{Float64,Nothing}
    err_float128::Union{Float64,Nothing}
    seconds::Float64
    problems::Vector{String}
end

"""
    check_cell(cell; precision = 192, dir = datadir()) -> CellResult

Every check the replication makes on one double-precision rule:

1. the file is the deposited file (sha256 against `data/<family>/SHA256SUMS`);
2. it holds `N` nodes in `d` dimensions, `N` and `d` being those of the file name and header;
3. all weights are positive, they sum to 1 (to `1e-14`), and for the uniform weight every node
   lies in the open cube;
4. the rule is exact for **all** monomials of total degree at most `p`: the largest relative
   error, in `precision`-bit arithmetic on the double-precision numbers, is below the paper's
   gate `1e-11`, and equals the error the deposit's header claims;
5. the row of the paper's table — `N`, `ρ`, the error, Möller's bound, the shading — is
   reproduced character for character.  The error in that column is the one of the rule's
   extended-precision file, taken here from the deposit's record ([`recorded_extended`](@ref));
   `replicate(extended = …)` recomputes it from the file itself.
"""
function check_cell(c::Cell; precision::Int = 192, dir::AbstractString = datadir(),
                    rows = paper_rows(c.family; dir), sums = checksums(c.family; dir),
                    recorded = recorded_extended(c.family; dir))
    t0 = time(); problems = String[]
    rel = "rules/" * basename(c.file)
    haskey(sums, rel) ? (sums[rel] == sha256file(c.file) || push!(problems, "sha256 differs from the deposit's SHA256SUMS")) :
                        push!(problems, "no checksum on record")
    X, w = load_rule64(c.file; precision)
    size(X) == (c.n, c.d) || push!(problems, "file holds $(size(X, 1)) nodes in $(size(X, 2)) dimensions, not $(c.n) in $(c.d)")
    v = verify_rule(X, w, c.p, FAMILY[c.family]; precision)
    v.minw > 0 || push!(problems, "a weight is not positive (min $(v.minw))")
    v.sumwdev ≤ 1e-14 || push!(problems, @sprintf("weights do not sum to 1 (off by %.2e)", v.sumwdev))
    c.family == "le" && !(v.minnode > 0 && v.maxnode < 1) && push!(problems, "a node lies outside the open cube")
    v.err ≤ GATE || push!(problems, @sprintf("relative error %.3e exceeds the gate %.0e", v.err, GATE))
    claimed = header_error(c.file)
    claimed === nothing || isapprox(claimed, v.err; rtol = 1e-5, atol = 1e-300) ||
        push!(problems, @sprintf("verified error %.6e, the file's header claims %.6e", v.err, claimed))
    pr = get(rows, (c.d, c.p), nothing)
    terr = get(recorded, (c.d, c.p, c.n), nothing)
    terr === nothing && push!(problems, "the deposit records no extended-precision error for this cell; the table row is built from the double-precision error")
    row = table_row(c, something(terr, v.err), pr === nothing ? nothing : pr.prev, pr === nothing ? "---" : pr.src;
                    symfloor = pr === nothing || c.family ≠ "gh" ? nothing : pr.pair)
    paper = pr === nothing ? nothing : paper_row_tex(pr)
    if pr === nothing
        push!(problems, "cell is not in the paper's table")
    else
        paper == row || push!(problems, "table row differs: computed `$row`, paper `$paper`")
        pr.prev !== nothing && pr.shade == "" && c.n ≥ pr.prev && push!(problems, "shown as an improvement, but N ≥ prev")
    end
    CellResult(c, v.err, v.minw, v.sumwdev, row, paper, nothing, nothing, time() - t0, problems)
end

"""
    check_extended(cell, extfile; float64file = cell.file) -> (err, problems)

Checks on an extended-precision file of the deposits (`rules_extended/*.mp.csv`): weights
positive, nodes inside the cube for the uniform weight, exact to the target of its family
(`1e-68` for both weights) in arithmetic wide enough for its digits, and — the link
between the two files of a cell — rounding it to double precision gives the double-precision
file, row for row and bit for bit.
"""
function check_extended(c::Cell, extfile::AbstractString; float64file::AbstractString = c.file)
    problems = String[]
    m = match(r"\((\d+) significant digits\)", something(header(extfile, "precision shipped"), ""))
    digits = m === nothing ? 40 : parse(Int, m[1])
    prec = max(256, ceil(Int, (digits + 45) * 3.33))
    X, w = load_rule(extfile; precision = prec)
    v = verify_rule(X, w, c.p, FAMILY[c.family]; precision = prec)
    v.minw > 0 || push!(problems, "extended file: a weight is not positive")
    c.family == "le" && !(v.minnode > 0 && v.maxnode < 1) && push!(problems, "extended file: a node lies outside the open cube")
    v.err ≤ EXT_TARGET[c.family] || push!(problems, @sprintf("extended file: error %.3e misses the target %.0e", v.err, EXT_TARGET[c.family]))
    F = [parse.(Float64, split(l, ',')) for l in eachline(float64file) if isdata(l)]
    R = [vcat(Float64.(X[i, :]), Float64(w[i])) for i in 1:size(X, 1)]
    F == R || push!(problems, "extended file does not round to the double-precision file row for row")
    v.err, problems
end

"""
    check_float128(cell, extfile) -> (err, problems)

The extended-precision file rounded to IEEE binary128 (`Float128`: 113-bit significand, machine
epsilon `2^-112 ≈ 1.93e-34`).  Every number is read from its decimal string straight into a
113-bit `BigFloat`, which MPFR rounds correctly, so the values are bit for bit those a program
holds after parsing the file into `Float128`; no quadruple-precision package is needed.  The
rounded rule must have positive weights and be exact over all monomials of degree at most `p` to
`10` machine epsilons — a few epsilons is the floor of the format, as it is for the
double-precision files — and the error must equal the one the file's header claims, where
the header carries that line (deposits built from 2026-09-19 on).
"""
function check_float128(c::Cell, extfile::AbstractString)
    problems = String[]
    m = match(r"\((\d+) significant digits\)", something(header(extfile, "precision shipped"), ""))
    prec = max(256, ceil(Int, ((m ≡ nothing ? 40 : parse(Int, m[1])) + 45) * 3.33))
    X, w = load_rule(extfile; precision = 113)
    X = BigFloat[BigFloat(x; precision = prec) for x ∈ X]; w = BigFloat[BigFloat(x; precision = prec) for x ∈ w]
    v = verify_rule(X, w, c.p, FAMILY[c.family]; precision = prec)
    v.minw > 0 || push!(problems, "Float128 rounding: a weight is not positive")
    v.err ≤ F128_BOUND || push!(problems, @sprintf("Float128 rounding: error %.3e exceeds 10 machine epsilons (%.2e)", v.err, F128_BOUND))
    claimed = nothing
    for l ∈ eachline(extfile)
        startswith(l, "#") || break
        h = match(r"^# rounded to IEEE binary128[^:]*: verified max relative monomial error (\S+)", l)
        h ≡ nothing || (claimed = parse(Float64, h[1]))
    end
    claimed ≡ nothing || isapprox(claimed, v.err; rtol = 1e-5, atol = 1e-300) ||
        push!(problems, @sprintf("Float128 rounding: verified error %.6e, the file's header claims %.6e", v.err, claimed))
    v.err, problems
end

"""
    replicate(; families = ("gh", "le"), maxnodes = typemax(Int), extended = nothing,
                outdir = "replication_output", precision = 192, verbose = true) -> Vector{CellResult}

Replicate the tables of the paper from the rules.

For every cell of the families asked for (and at most `maxnodes` nodes) this runs `check_cell`
— integrity, positivity, mass, exactness over all monomials in `precision`-bit arithmetic, and
a character-for-character rebuild of the cell's row in the paper's table — and writes to
`outdir`

* `replication_report.md`: the verdict, the headline counts and every discrepancy;
* `cells.csv`: one line per cell with the computed quantities;
* `best_known_gh_rows.tex`, `best_known_le_rows.tex`: the rebuilt table rows.

`extended` may name a directory that holds the two Zenodo deposits unpacked side by side
(`<extended>/gh/rules_extended`, `<extended>/le/rules_extended`); each extended-precision file
found there is then checked with `check_extended` and, rounded to quadruple precision, with
`check_float128`.

The full run takes minutes (most of it in the three largest `d = 5` cells) and uses all the
threads Julia was started with.  `maxnodes = 500` is a quick pass over the small cells.
"""
function replicate(; families = ("gh", "le"), maxnodes::Int = typemax(Int), extended = nothing,
                   outdir::AbstractString = "replication_output", precision::Int = 192,
                   verbose::Bool = true, dir::AbstractString = datadir())
    results = CellResult[]
    for fam in families
        rows = paper_rows(fam; dir); sums = checksums(fam; dir); recorded = recorded_extended(fam; dir)
        ext = Dict{Tuple{Int,Int},Cell}()
        extended === nothing || for e in cells(fam; dir = extended, sub = "rules_extended"); ext[(e.d, e.p)] = e; end
        for c in cells(fam; dir)
            c.n ≤ maxnodes || continue
            r = check_cell(c; precision, dir, rows, sums, recorded)
            if haskey(ext, (c.d, c.p))
                e = ext[(c.d, c.p)]
                eerr, eprob = e.n == c.n ? check_extended(c, e.file) : (nothing, ["extended file is for N = $(e.n)"])
                qerr, qprob = e.n == c.n ? check_float128(c, e.file) : (nothing, String[])
                # the table's error column is this file's error: what the deposit records must be what the file measures
                rec = get(recorded, (c.d, c.p, c.n), nothing)
                eerr === nothing || rec === nothing || isapprox(eerr, rec; rtol = 1e-5, atol = 1e-300) ||
                    push!(eprob, @sprintf("extended file measures %.6e, the deposit records %.6e (the paper's table prints the latter)", eerr, rec))
                r = CellResult(r.cell, r.err, r.minw, r.sumwdev, r.row, r.paper, eerr, qerr, r.seconds, vcat(r.problems, eprob, qprob))
            end
            push!(results, r)
            verbose && println(rpad(string(c), 24), @sprintf("err %.1e", r.err),
                               r.err_extended === nothing ? "" : @sprintf("  extended %.1e", r.err_extended),
                               r.err_float128 ≡ nothing ? "" : @sprintf("  Float128 %.1e", r.err_float128),
                               isempty(r.problems) ? "  ok" : "  PROBLEM: " * join(r.problems, "; "),
                               @sprintf("  (%.1f s)", r.seconds))
        end
        # a cell of the paper's table with no rule behind it is a failure too
        have = Set((r.cell.d, r.cell.p) for r in results if r.cell.family == fam)
        maxnodes == typemax(Int) && for k in sort(collect(keys(rows)))
            k ∈ have || push!(results, CellResult(Cell(fam, k[1], k[2], rows[k].N, ""), NaN, NaN, NaN, "", paper_row_tex(rows[k]), nothing, nothing, 0.0,
                                                  ["the paper's table has this cell, the data have no rule for it"]))
        end
    end
    write_report(results, outdir; families, maxnodes, extended, precision)
    verbose && println(all(r -> isempty(r.problems), results) ? "REPLICATED: " : "NOT REPLICATED: ",
                       length(results), " cells, ", count(r -> !isempty(r.problems), results), " with problems; report in ", outdir)
    results
end

"True when every cell passed every check."
replicated(results::Vector{CellResult}) = !isempty(results) && all(r -> isempty(r.problems), results)

function write_report(results, outdir; families, maxnodes, extended, precision)
    mkpath(outdir)
    open(joinpath(outdir, "cells.csv"), "w") do io
        println(io, "family,d,p,N,rho,moller_bound,rel_err_float64,rel_err_extended,rel_err_float128,min_weight,sum_weights_minus_1,seconds,problems")
        for r in results
            c = r.cell
            println(io, join((c.family, c.d, c.p, c.n, @sprintf("%.6f", rho(c.n, c.d, c.p)), moller_bound(c.d, c.p),
                              @sprintf("%.6e", r.err), r.err_extended === nothing ? "" : @sprintf("%.6e", r.err_extended),
                              r.err_float128 ≡ nothing ? "" : @sprintf("%.6e", r.err_float128),
                              @sprintf("%.6e", r.minw), @sprintf("%.3e", r.sumwdev), @sprintf("%.1f", r.seconds),
                              "\"" * replace(join(r.problems, "; "), "\"" => "'") * "\""), ","))
        end
    end
    for fam in families
        open(joinpath(outdir, "best_known_$(fam)_rows.tex"), "w") do io
            d = 0
            for r in results
                r.cell.family == fam || continue
                r.cell.d == d || (d = r.cell.d; println(io, "\\multicolumn{7}{@{}l}{\$d = $d\$} \\\\"))
                println(io, r.row)
            end
        end
    end
    open(joinpath(outdir, "replication_report.md"), "w") do io
        bad = [r for r in results if !isempty(r.problems)]
        println(io, "# Replication report\n")
        println(io, isempty(bad) ? "**Replicated.**" : "**Not replicated: $(length(bad)) of $(length(results)) cells have problems.**",
                "  $(length(results)) cells checked in $(precision)-bit arithmetic",
                maxnodes == typemax(Int) ? "" : " (cells with at most $maxnodes nodes only)",
                extended === nothing ? "; extended-precision files not checked." : "; extended-precision files checked for $(count(r -> r.err_extended !== nothing, results)) cells.", "\n")
        println(io, "| family | cells | at Möller's bound | below the literature | equal to the literature | no earlier rule besides the product grid | largest error |")
        println(io, "|---|---|---|---|---|---|---|")
        for fam in families
            R = [r for r in results if r.cell.family == fam && r.row ≠ ""]
            rows = paper_rows(fam)
            prev(r) = (pr = get(rows, (r.cell.d, r.cell.p), nothing); pr === nothing ? nothing : pr.prev)
            println(io, "| ", FAMILY_NAME[fam], " | ", length(R), " | ", count(r -> r.cell.n == moller_bound(r.cell.d, r.cell.p), R),
                    " | ", count(r -> prev(r) !== nothing && r.cell.n < prev(r), R),
                    " | ", count(r -> prev(r) !== nothing && r.cell.n == prev(r), R),
                    " | ", count(r -> prev(r) === nothing, R),
                    " | ", isempty(R) ? "" : @sprintf("%.1e", maximum(r.err for r in R)), " |")
        end
        if !isempty(bad)
            println(io, "\n## Problems\n")
            for r in bad; println(io, "- **", r.cell, "**: ", join(r.problems, "; ")); end
        end
        println(io, "\nChecks per cell: file checksum; N and d; positive weights summing to 1; nodes inside the cube (uniform weight); ",
                "relative error over all monomials of total degree ≤ p below 1e-11 and equal to the deposit's claim; ",
                "the paper's table row rebuilt character for character.")
        q = [r.err_float128 for r ∈ results if r.err_float128 ≢ nothing]
        isempty(q) || println(io, @sprintf("\nExtended-precision files rounded to IEEE binary128 (Float128): largest error %.2e, that is %.2f machine epsilons (2^-112); the bound checked is 10.", maximum(q), maximum(q) / EPS128))
    end
end
