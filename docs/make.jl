using Documenter
using Modulo2

DocMeta.setdocmeta!(Modulo2, :DocTestSetup, quote
        using Modulo2
    # for jldoctest in docstrings
    end; recursive = true)

makedocs(sitename = "Modulo2.jl",
    modules = [Modulo2],
    format = Documenter.HTML(),
    warnonly = true)
