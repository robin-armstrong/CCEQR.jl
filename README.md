This package provides `cceqr!`, an Julia implementation of the "collect, commit, expand" column-pivoted QR factorization. See Armstrong and Damle, *Collect, Commit, Expand: Efficient CPQR-Based Column Selection for Extremely Wide Matrices*, 2024. To install this package, run the following commands in a Julia terminal:
```
import Pkg
Pkg.install("https://github.com/robin-armstrong/CCEQR.git")
```
After this, run `using CCEQR` to bring `cceqr!` into the namespace. Running `?cceqr!` will show detailed documentation for the function.