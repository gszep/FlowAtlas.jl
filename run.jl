using Pkg
Pkg.activate(".")

if isfile("cytometry.so")
    using Pluto
    Pluto.run(notebook="cytometry.jl",sysimage="cytometry.so")
else
    println("Running app for the first time requires building cytometry.so\nThis may take minutes, feel free to grab a cuppa tea ;)")
    include("build/sysimage.jl")
end