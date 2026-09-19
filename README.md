# QuadricepsReplicationPackage.jl

Replication package for

> J. Pinkse, *Positive weight Hermite and Legendre quadrature rules*, 2026.

The paper's contribution is two tables of positive-weight quadrature rules of odd degree of
exactness `p`, in `d = 2, …, 5` dimensions: 57 rules for the Gaussian weight `N(0, I_d)` (GH)
and 85 for the uniform weight on `[0,1]^d` (Le). This package holds every one of those rules in
double precision and rebuilds the tables from them.

```julia
using QuadricepsReplicationPackage

results = replicate()      # all 142 cells
replicated(results)        # true when every check passed
```

Start Julia with threads (`julia -t 8`): the work is exact arithmetic on up to 13,199 nodes and
65,780 monomials per rule, and the full run takes about seven minutes on 16 threads, extended-precision files included. `replicate(maxnodes = 500)` is a
quick pass over the small cells.

## What is replicated

For each cell `replicate` checks that

1. the rule file is the deposited file (sha256 against the deposit's `SHA256SUMS`);
2. it has the `N` nodes and `d` dimensions the tables say;
3. all weights are positive and sum to 1, and for the uniform weight every node lies in the
   open cube;
4. the rule is exact for **all** monomials of total degree at most `p`: the largest relative
   error

       max_{|a| ≤ p} |Σ_s w_s x_s^a − E x^a| / max(Σ_s |w_s| |x_s^a|, 1),

   evaluated in 192-bit arithmetic on the double-precision numbers, is below the paper's
   acceptance threshold `1e-11` and equals the figure in the file's header;
5. the cell's row in the paper's table — `N`, `ρ = N^{1/d}/q` with `q = (p+1)/2`, the error,
   Möller's lower bound, the shading — is reproduced character for character
   (`data/paper/best_known_gh.tex`, `best_known_le.tex` are the tables as typeset).

It writes `replication_report.md`, `cells.csv` and the rebuilt table rows to
`replication_output/`.

Two columns of the tables describe the literature, not the rules: `prev`, the best node count
published before the paper, and `src`, where it was published. They are taken from the paper's
tables as given; the sources are listed in the paper's bibliography.

### The extended-precision files

The Zenodo deposits also hold every rule to 40 digits (GH) or 80 digits (Le). Those files are
not in this package (32 MB). With the two deposits unpacked side by side,

```julia
replicate(extended = "/path/to/deposits")   # expects gh/rules_extended and le/rules_extended there
```

additionally checks, for each such file, positivity, exactness to `1e-34` (GH) or `1e-68` (Le)
in arithmetic wide enough for its digits, and that rounding it to double precision gives the
double-precision file row for row and bit for bit.

## What is not replicated

The search that found the rules. It ran for weeks on two clusters and several desktops, from
random starts and under wall-clock budgets, and is not deterministic. Nothing in the tables
depends on it: a rule proves its own node count, exactness and positivity, which is what the
checks above establish, and the only minimality statements in the paper are the cells that
attain Möller's bound, which `moller_bound` computes. The solvers are in a separate
repository, PositiveWeightQuadratureSolvers.jl; the rules themselves, for use rather than for
checking, are in Quadriceps.jl.

## Functions

| | |
|---|---|
| `replicate(; families, maxnodes, extended, outdir, precision)` | run every check and write the report |
| `replicated(results)` | `true` when every cell passed |
| `check_cell(cell)` | the checks on one double-precision rule |
| `check_extended(cell, file)` | the checks on one extended-precision file |
| `verify_rule(nodes, weights, p, family; precision)` | largest relative monomial error, weight and node diagnostics |
| `moller_bound(d, p)` | Möller's lower bound on the number of nodes |
| `rho(N, d, p)` | the tables' quality measure |
| `cells(family)`, `load_rule64(path)`, `load_rule(path)` | the rule files and how to read them |
| `paper_rows(family)`, `table_row(cell, err, prev, src)` | the paper's rows, and a rebuilt row |

## Data

`data/<gh|le>/rules/` are the double-precision files of the Zenodo deposits, unchanged: comment
lines with the cell, the verified error, the credit and its evidence, then one node per line
(`x1,…,xd,w`). `summary.csv` is the deposit's catalog and `SHA256SUMS` its checksums for these
files. GH nodes are in the frame of the standard normal density; Le nodes are in `[0,1]^d`;
weights sum to 1 in both.

Some rules are, or descend from, rules published by others; see [`NOTICE.md`](NOTICE.md) and
the `# credit:` line of each file.

## License

The code is under the MIT license ([`LICENSE`](LICENSE)). The rules derived from Diallo and
Worku's data carry their MIT notice, reproduced in `NOTICE.md`.
