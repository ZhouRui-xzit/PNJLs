using Pkg


needs = [
    
    "ForwardDiff", #computational_MFE 
    "NLsolve",
    "FastGaussQuadrature",
    "Peaks",
    "QuadGK",
    "DataInterpolations",
    "Roots",
    "FiniteDifferences",
    "NaNMath",
    "SpecialFunctions",
    "ApproxFun",# pdes
    "LinearAlgebra",
    "ModelingToolkit",
    "DifferentialEquations",
    "MethodOfLines",
    "MeshGrid",

    "DataFrames",# data and plot
    "CSV",
    "Plots",
    "LaTeXStrings",
    
    "BenchmarkTools",# others
    "Revise",

    "LanguageServer", # for terminal editor like neovim
    "SymbolServer",
    "StaticLint",
    "CSTParser",
]

for pkg in needs
    println("$(pkg) is installing...")
    Pkg.add(pkg)
end