using QuadricepsReplicationPackage
using Test

const QRP = QuadricepsReplicationPackage

@testset "QuadricepsReplicationPackage" begin
    @testset "Möller's bound" begin
        # closed forms in Orive, Santos-León and Spalević (2020): 2d at degree 3, d² + d + 1 at degree 5
        for d in 2:8
            @test moller_bound(d, 3) == 2d
            @test moller_bound(d, 5) == d^2 + d + 1
        end
        @test moller_bound(2, 1) == 1
        @test moller_bound(2, 11) == 24
        @test_throws ArgumentError moller_bound(3, 4)
        # every bound printed in the paper's tables
        for fam in ("gh", "le"), (k, r) in paper_rows(fam)
            @test r.moller == moller_bound(k...)
        end
    end

    @testset "the pair floor" begin
        # 2⌈M/(d+1)⌉ by hand at d = 3, p = 21: M = 946 even-degree monomials, 4 unknowns per pair
        @test pair_floor(3, 21) == 474
        @test pair_floor(5, 21) == 9820
        @test pair_floor(2, 1) == 1
        @test_throws ArgumentError pair_floor(3, 4)
        # the pair column is the uniform-weight table's, and no rule there is below it (the Gaussian
        # rules carry more symmetry and several are, which is why that table prints sym instead)
        for (k, r) in paper_rows("le")
            @test r.pair == pair_floor(k...)
            @test r.N ≥ r.pair
        end
        # the Gaussian table's floor column is sym, taken as given; the paper's claim about it is
        # that no rule has fewer nodes than sym, and sym is never below Möller's bound
        for (_, r) in paper_rows("gh")
            @test r.pair !== nothing && r.N ≥ r.pair ≥ r.moller
        end
    end

    @testset "verification kernel" begin
        # the 2-point Gauss rules, exact to degree 3, in one and two dimensions
        b(x) = BigFloat(x; precision = 192)
        g = 1 / sqrt(b(3))
        X1 = reshape([(1 - g) / 2, (1 + g) / 2], 2, 1); w1 = [b(1) / 2, b(1) / 2]
        @test verify_rule(X1, w1, 3, "legendre").err < 1e-55
        @test verify_rule(X1, w1, 5, "legendre").err > 1e-3            # not exact to degree 5
        X2 = BigFloat[x for x in (-1, 1), _ in 1:1]
        @test verify_rule(X2, w1, 3, "hermite").err < 1e-55
        Xt = BigFloat[s * b(1) for (s, t) in ((-1, -1), (-1, 1), (1, -1), (1, 1)), _ in 1:1]
        Xt = hcat(Xt, BigFloat[t * b(1) for (s, t) in ((-1, -1), (-1, 1), (1, -1), (1, 1))])
        v = verify_rule(Xt, fill(b(1) / 4, 4), 3, "hermite")
        @test v.err < 1e-55 && v.minw == 0.25 && v.sumwdev == 0
        @test_throws ArgumentError verify_rule(X1, w1, 3, "laguerre")
    end

    @testset "data" begin
        for (fam, ncell) in (("gh", 57), ("le", 85))
            cs = cells(fam)
            @test length(cs) == ncell
            @test Set((c.d, c.p) for c in cs) == Set(keys(paper_rows(fam)))
            @test all(isodd(c.p) && 2 ≤ c.d ≤ 5 for c in cs)
        end
    end

    @testset "small cells replicate" begin
        out = mktempdir()
        res = replicate(; maxnodes = 120, outdir = out, verbose = false)
        @test length(res) ≥ 40
        for r in res
            @test isempty(r.problems)
        end
        @test replicated(res)
        @test isfile(joinpath(out, "replication_report.md"))
        @test occursin("**Replicated.**", read(joinpath(out, "replication_report.md"), String))
        @test countlines(joinpath(out, "cells.csv")) == length(res) + 1
    end

    @testset "a tampered rule is caught" begin
        c = first(c for c in cells("gh") if c.d == 2 && c.p == 5)
        tmp = mktempdir(); dst = joinpath(tmp, "gh", "rules"); mkpath(dst)
        lines = readlines(c.file)
        i = findlast(QRP.isdata, lines)
        f = split(lines[i], ','); f[end] = string(parse(Float64, f[end]) * (1 + 1e-9)); lines[i] = join(f, ',')
        bad = joinpath(dst, basename(c.file)); write(bad, join(lines, "\n") * "\n")
        r = check_cell(Cell(c.family, c.d, c.p, c.n, bad))
        @test any(occursin("sha256", p) for p in r.problems)
        @test any(occursin("header claims", p) || occursin("sum to 1", p) for p in r.problems)
    end

    @testset "extended-precision files" begin
        ext = joinpath(@__DIR__, "data")
        for e in cells(; dir = ext, sub = "rules_extended")
            c = first(c for c in cells(e.family) if (c.d, c.p, c.n) == (e.d, e.p, e.n))
            err, problems = check_extended(c, e.file)
            @test isempty(problems)
            @test err ≤ 1e-68                                    # both weights ship 80 digits (GH since 2026-09-20)
            rec = QRP.recorded_extended(e.family)[(e.d, e.p, e.n)]  # the figure the paper's table prints
            @test isapprox(err, rec; rtol = 1e-5)
            err128, problems128 = check_float128(c, e.file)     # the same file rounded to IEEE binary128
            @test isempty(problems128)
            @test err < err128 ≤ 10 * 2.0^-112
            @test occursin("rounded to IEEE binary128", read(e.file, String))   # the deposit states the figure; check_float128 compared it
        end
        @test length(cells(; dir = ext, sub = "rules_extended")) == 2
    end
end
