# from https://github.com/JuliaLang/julia/pull/36425
function detectwsl()
    Sys.islinux() &&
    isfile("/proc/sys/kernel/osrelease") &&
    occursin(r"Microsoft|WSL"i, read("/proc/sys/kernel/osrelease", String))
end

function browser(url::AbstractString)
    if Sys.isapple()
        Base.run(`open $url`)

    elseif Sys.iswindows() || detectwsl()
        Base.run(`powershell.exe Start $url`)

    elseif Sys.islinux()
        Base.run(`xdg-open $url`)

    else
        printstyled("[FlowAtlas] ", color = :blue)
        println("Could not open default browser")
    end
end