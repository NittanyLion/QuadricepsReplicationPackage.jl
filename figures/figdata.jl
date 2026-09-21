# figdata.jl — data for the two figures of quadriceps.tex (fig_nodes.tex, fig_rho.tex).
#
#   julia figdata.jl [<output folder>]        (Julia 1.13; a second or two)
#
# The same file runs in two places.  In the author's project folder (paper_tables/) it reads the
# deposits under ../publish and the two generated tables next to itself, and writes fig/*.dat next
# to itself.  Inside QuadricepsReplicationPackage.jl (figures/) it reads the package's data/, writes
# to replication_output/figures/, and ends by comparing what it wrote with the reference output
# shipped in figures/fig/, which is what the paper's figures were drawn from.
#
#   fig/rho_<gh|le>_d<d>.dat   ρ = N^(1/d)/q against p, with the counting floor and Möller's bound on
#                              the same scale, read off best_known_<gh|le>.tex
#   fig/nodes_<gh|le>.dat      nodes and weights of the largest planar rule of each table
#
# No package beyond the standard library is needed.
using Printf, DelimitedFiles

const HERE   = @__DIR__
const ROOT   = dirname(HERE)
const PACKAGED = isdir(joinpath(ROOT, "data", "gh", "rules"))       # inside the replication package
rules_dir(sub) = PACKAGED ? joinpath(ROOT, "data", sub, "rules") : joinpath(ROOT, "publish", sub, "rules")
const TABLES = PACKAGED ? joinpath(ROOT, "data", "paper") : HERE      # best_known_<gh|le>.tex
const OUT    = let a = filter(!startswith("--"), ARGS)
    !isempty(a) ? abspath(a[1]) : PACKAGED ? joinpath(ROOT, "replication_output", "figures") : HERE
end
const FIG    = joinpath(OUT, "fig")
mkpath(FIG)

# ρ against p, with the counting floor and Möller's bound on the same scale, read off the paper's tables
for (fam, file) ∈ (("gh", "best_known_gh.tex"), ("le", "best_known_le.tex"))
    d = 0; io = nothing; seen = Set{Int}()
    for l ∈ eachline(joinpath(TABLES, file))
        m = match(r"^\\multicolumn.*\$d = (\d)\$", l)      # block headers, "(cont.)" ones included
        if m ≢ nothing
            io ≡ nothing || close(io)
            d = parse(Int, m[1]); io = open(joinpath(FIG, "rho_$(fam)_d$(d).dat"), d ∈ seen ? "a" : "w")
            d ∈ seen || println(io, "p rho floor moller"); push!(seen, d)
            continue
        end
        c = strip.(split(replace(l, r"\\rowcolor\{[^}]*\}" => "", "\\\\" => ""), "&"))
        (length(c) == 8 && all(isdigit, c[1])) || continue
        p, n, mo, fl = parse.(Int, (c[1], c[2], c[5], c[6])); q = (p + 1) ÷ 2
        @printf(io, "%d %.5f %.5f %.5f\n", p, n^(1 / d) / q, fl^(1 / d) / q, mo^(1 / d) / q)
    end
    io ≡ nothing || close(io)
end
# two planar rules, nodes and weights
for (fam, sub, pat) ∈ (("gh", "gh", "hermite_d2_p33_"), ("le", "le", "legendre_d2_p77_"))
    dir = rules_dir(sub); f = only(filter(startswith(pat), readdir(dir)))
    A = readdlm(joinpath(dir, f), ',', Float64; comments = true, comment_char = '#', header = true)[1]
    open(joinpath(FIG, "nodes_$(fam).dat"), "w") do io
        println(io, "x y w logw")
        for i ∈ axes(A, 1); @printf(io, "%.8f %.8f %.6e %.4f\n", A[i, 1], A[i, 2], A[i, 3], log10(A[i, 3])); end
    end
end

# Comparison with the reference, only when the output went somewhere else than next to this file.
if OUT ≠ HERE && isdir(joinpath(HERE, "fig"))
    bad = String[]; n = 0
    for f ∈ readdir(joinpath(HERE, "fig"))
        endswith(f, ".dat") || continue
        g = joinpath(FIG, f); isfile(g) || (push!(bad, "$f: not produced"); continue)
        A, B = readdlm(joinpath(HERE, "fig", f); skipstart = 1), readdlm(g; skipstart = 1)
        size(A) == size(B) || (push!(bad, "$f: $(size(B, 1)) rows, reference has $(size(A, 1))"); continue)
        all(abs(a - b) ≤ 1e-3 * max(abs(a), abs(b)) for (a, b) ∈ zip(A, B)) || push!(bad, "$f: values differ")
        global n += 1
    end
    println(isempty(bad) ? "REPLICATED: all $n figure data files agree with the reference in figures/fig" :
            "DIFFERENT from the reference ($(length(bad)) of $n files):\n  " * join(bad, "\n  "))
    exit(isempty(bad) ? 0 : 1)
end
