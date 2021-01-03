using Pkg
Pkg.activate(".")

printstyled(color=:cyan,"[pkg] ")
printstyled("Installing/updating julia package dependencies...\n")

Pkg.instantiate()
Pkg.update()

if ~isfile("cytometry.so")
    printstyled(color=:yellow,"[sysimage] ")
    printstyled("Building cytometry.so\nThis may take up to 30mins, feel free to grab a cuppa tea ;)\n")
    include("build/sysimage.jl")
end

printstyled(color=:blue,"[Pluto] ")
printstyled("Launching notebook...\n")
using Pluto
Pluto.run(notebook="cytometry.jl",sysimage="cytometry.so")