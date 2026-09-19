const FAMILY = Dict("gh" => "hermite", "le" => "legendre")
const FAMILY_NAME = Dict("gh" => "Gaussian weight N(0, I_d)", "le" => "uniform weight on [0,1]^d")

"The package's data directory (the double-precision rules, the paper's tables, the checksums)."
datadir(parts...) = normpath(joinpath(@__DIR__, "..", "data", parts...))

"""
    Cell

One cell of the paper's tables: `family` is `"gh"` or `"le"`, `d` the dimension, `p` the degree
of exactness, `n` the number of nodes, `file` the rule file.
"""
struct Cell
    family::String
    d::Int
    p::Int
    n::Int
    file::String
end
Base.show(io::IO, c::Cell) = print(io, uppercasefirst(c.family), " d=", c.d, " p=", c.p, " N=", c.n)

const RULEFILE = r"^(hermite|legendre)_d(\d+)_p(\d+)_q(\d+)_n(\d+)\.(?:mp\.)?csv$"

"""
    cells(family = nothing; dir = datadir()) -> Vector{Cell}

The rule files under `dir/<family>/rules`, sorted by dimension and degree.  `family` is `"gh"`,
`"le"` or `nothing` for both.
"""
function cells(family = nothing; dir::AbstractString = datadir(), sub::AbstractString = "rules")
    out = Cell[]
    for fam in (family === nothing ? ("gh", "le") : (family,))
        root = joinpath(dir, fam, sub)
        isdir(root) || continue
        for f in readdir(root)
            m = match(RULEFILE, f)
            m === nothing && continue
            m[1] == FAMILY[fam] || continue
            push!(out, Cell(fam, parse(Int, m[2]), parse(Int, m[3]), parse(Int, m[5]), joinpath(root, f)))
        end
    end
    sort!(out; by = c -> (c.family, c.d, c.p))
end

isdata(l) = occursin(r"^\s*[-+0-9.]", l)

"""
    load_rule(path; precision = 192) -> (nodes, weights)

Read a rule file: comment lines start with `#`, one header line names the columns, and every
other line holds the `d` coordinates of a node followed by its weight.  Numbers are read as
decimal strings straight into `BigFloat` at `precision` bits.  For a double-precision file that
is the exact value of the `Float64` whenever the file holds its shortest round-trip
representation — use `load_rule64` to go through `Float64` explicitly.
"""
function load_rule(path::AbstractString; precision::Int = 192)
    rows = [[BigFloat(String(s); precision = precision) for s in split(l, ',')] for l in eachline(path) if isdata(l)]
    isempty(rows) && error("no data lines in $path")
    d = length(rows[1]) - 1
    X = BigFloat[rows[i][k] for i in eachindex(rows), k in 1:d]
    w = BigFloat[rows[i][d + 1] for i in eachindex(rows)]
    X, w
end

"""
    load_rule64(path; precision = 192) -> (nodes, weights)

Read a double-precision rule file the way a user of the rule would: every number is parsed as a
`Float64` and that `Float64` is converted exactly to `BigFloat`.  The verified error is then the
error of the rule as a program holds it, with no help from digits beyond double precision.
"""
function load_rule64(path::AbstractString; precision::Int = 192)
    rows = [parse.(Float64, split(l, ',')) for l in eachline(path) if isdata(l)]
    isempty(rows) && error("no data lines in $path")
    d = length(rows[1]) - 1
    X = BigFloat[BigFloat(rows[i][k]; precision = precision) for i in eachindex(rows), k in 1:d]
    w = BigFloat[BigFloat(rows[i][d + 1]; precision = precision) for i in eachindex(rows)]
    X, w
end

"The value of a `# key: value` header line of a rule file, or `nothing`."
function header(path::AbstractString, key::AbstractString)
    for l in eachline(path)
        startswith(l, "#") || break
        m = match(Regex("^#\\s*" * key * ":\\s*(.*)\$"), l)
        m === nothing || return String(strip(m[1]))
    end
    nothing
end

"The error the deposit's own header claims for a rule file, or `nothing`."
function header_error(path::AbstractString)
    for l in eachline(path)
        startswith(l, "#") || break
        m = match(r"^# verified max relative monomial error[^:]*:\s*(\S+)", l)
        m === nothing || return parse(Float64, m[1])
    end
    nothing
end

sha256file(path) = bytes2hex(open(sha256, path))

"The checksums shipped with the data: `relative path => sha256`, for family `fam`."
function checksums(fam::AbstractString; dir::AbstractString = datadir())
    out = Dict{String,String}()
    f = joinpath(dir, fam, "SHA256SUMS")
    isfile(f) || return out
    for l in eachline(f)
        m = match(r"^([0-9a-f]{64})\s+\*?(\S.*)$", l)
        m === nothing || (out[String(m[2])] = String(m[1]))
    end
    out
end
