using Documenter
using Incertus

makedocs(
    sitename = "Incertus",
    format = Documenter.HTML(),
    modules = [Incertus]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/yevgenryeznik/Incertus.jl.git"
)
