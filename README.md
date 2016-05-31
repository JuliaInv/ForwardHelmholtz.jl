[![Build Status](https://travis-ci.org/JuliaInv/ForwardHelmholtz.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/ForwardHelmholtz.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/ForwardHelmholtz.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/ForwardHelmholtz.jl?branch=master)

# ForwardHelmholtz.jl
A package for defining and solving the Helmholtz equation using the shifted Laplacian multigrid solver. 
This is done using the geometric MG version in Multigrid.jl.

# Requirements

This package is intended to use with julia versions 0.4.x.

This package is an add-on for [`jInv`](https://github.com/JuliaInv/jInv.jl), which needs to be installed. 

# Installation

```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/Multigrid.jl","Multigrid")
Pkg.clone("https://github.com/JuliaInv/ForwardHelmholtz.jl","ForwardHelmholtz")
Pkg.clone("https://github.com/JuliaInv/ParSpMatVec.jl","ParSpMatVec")
Pkg.build("ParSpMatVec");

Pkg.test("ForwardHelmholtz")
```

# Examples
See //examples.
