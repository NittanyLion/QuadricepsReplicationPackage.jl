# update_data.jl — refresh data/ from the project's sync folder (not part of the package API).
#
#   julia build/update_data.jl [<sync folder>] [--check]
#
# Copies, for each family, the double-precision rules of the planned Zenodo deposit
# (publish/<gh|le>/rules/*.csv) with the deposit's summary.csv and the `rules/` lines of its
# SHA256SUMS, the paper's two tables (paper_tables/best_known_<gh|le>.tex), and two small
# extended-precision files for the test suite, and the paper's two figures (paper_tables/figdata.jl, the
# figure sources and the figure data → figures/).  A rule file that is no longer in the deposit is
# removed.  With --check nothing is written: the exit status is 1 when data/ is out of date.
# The extended-precision files themselves are not copied (32 MB); `replicate(extended = …)`
# reads them from the unpacked deposits.
const SYNC = let a = filter(!startswith("--"), ARGS)
    isempty(a) ? joinpath(homedir(), "Dropbox", "oldDesignedQuadrature-sync") : a[1]
end
const CHECK = "--check" ∈ ARGS
const PKG = normpath(joinpath(@__DIR__, ".."))
stale = String[]
function put(src, dst)
    if !isfile(dst) || read(src) ≠ read(dst)
        push!(stale, relpath(dst, PKG))
        CHECK || (mkpath(dirname(dst)); cp(src, dst; force = true))
    end
end
for fam in ("gh", "le")
    dep = joinpath(SYNC, "publish", fam); out = joinpath(PKG, "data", fam)
    isdir(joinpath(dep, "rules")) || error("no deposit at $dep")
    keep = Set{String}()
    for f in readdir(joinpath(dep, "rules"))
        endswith(f, ".csv") || continue
        put(joinpath(dep, "rules", f), joinpath(out, "rules", f)); push!(keep, f)
    end
    isdir(joinpath(out, "rules")) && for f in readdir(joinpath(out, "rules"))
        f ∈ keep && continue
        push!(stale, "data/$fam/rules/$f (gone from the deposit)")
        CHECK || rm(joinpath(out, "rules", f))
    end
    put(joinpath(dep, "summary.csv"), joinpath(out, "summary.csv"))
    sums = join((l for l in eachline(joinpath(dep, "SHA256SUMS")) if occursin(r"\s\*?rules/", l)), "\n") * "\n"
    tmp = tempname(); write(tmp, sums); put(tmp, joinpath(out, "SHA256SUMS")); rm(tmp)
    put(joinpath(SYNC, "paper_tables", "best_known_$fam.tex"), joinpath(PKG, "data", "paper", "best_known_$fam.tex"))
end
for (fam, f) in (("gh", "hermite_d2_p5_q3_n7.mp.csv"), ("le", "legendre_d3_p5_q3_n13.mp.csv"))
    src = joinpath(SYNC, "publish", fam, "rules_extended", f)
    isfile(src) && put(src, joinpath(PKG, "test", "data", fam, "rules_extended", f))
end
# The paper's two figures: the script that writes their data, their pgfplots sources, and the data as
# the paper used it (figures/fig/*.dat is what the figures are drawn from and what a rerun is compared with).
let src = joinpath(SYNC, "paper_tables"), out = joinpath(PKG, "figures")
    for f in ("figdata.jl", "figstyle.tex", "fig_nodes.tex", "fig_rho.tex")
        put(joinpath(src, f), joinpath(out, f))
    end
    keep = Set(f for f in readdir(joinpath(src, "fig")) if endswith(f, ".dat"))
    for f in keep; put(joinpath(src, "fig", f), joinpath(out, "fig", f)); end
    isdir(joinpath(out, "fig")) && for f in readdir(joinpath(out, "fig"))
        f ∈ keep && continue
        push!(stale, "figures/fig/$f (gone from the paper)")
        CHECK || rm(joinpath(out, "fig", f))
    end
end
println(isempty(stale) ? "data/ is up to date" : (CHECK ? "OUT OF DATE: " : "updated: ") * string(length(stale)) * " files" *
        (length(stale) ≤ 12 ? " (" * join(stale, ", ") * ")" : ""))
exit(CHECK && !isempty(stale) ? 1 : 0)
