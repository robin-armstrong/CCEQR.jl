This package provides the Julia functions `cceqr` and `cceqr!`, which implement a variant of the column-pivoted QR algorithm for column subset selection on extremely wide matrices. To install this package, run the following commands in a Julia terminal:
```
import Pkg
Pkg.add(url="https://github.com/robin-armstrong/CCEQR.jl.git")
```
After this, run `using CCEQR` to bring `cceqr` and `cceqr!` into the namespace. Running `?cceqr` or `?cceqr!` will show detailed documentation for the functions.

For information about the underlying algorithm, see R. Armstrong and A. Damle, *Collect, Commit, Expand: Efficient CPQR-Based Column Selection for Extremely Wide Matrices*, 2025, https://export.arxiv.org/abs/2501.18035.