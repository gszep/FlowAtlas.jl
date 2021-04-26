import Pkg
Pkg.activate(".")

printstyled(color = :cyan,"[pkg] ")
printstyled("Installing/updating julia package dependencies...\n")

Pkg.instantiate()
Pkg.update()

printstyled("[FlowAtlas] ", color = :blue)
println("Compling Julia code... please be patient :) Will release pre-built sysimages in future")

using JSServe: Session, @js_str, file_server, response_404, Routes, App, DOM, Observable, on
using JSServe, HTTP, JSON, Images, FileIO, ImageIO
using GigaScatter, ColorSchemes

using Serialization: serialize,deserialize
using FlowWorkspace:inpolygon

using FlowWorkspace, StaticArrays, DataFrames, OrderedCollections, Glob
using StatsBase, GigaSOM, TSne
using StatsBase:normalize 

include("lib/colors.jl")
include("lib/selection.jl")
include("lib/embed.jl")

include("lib/tile.jl")
include("lib/gate.jl")
include("lib/count.jl")

include("lib/sys.jl")
include("lib/map.jl")
include("lib/sidebar.jl")

#################### server defaults
port = 3141
url = "http://localhost:$port"

#################### data inputs
workspace = "data/workspace.wsp"
files = glob"data/*/*.cleaned.fcs"
channelMap = Dict([

    "FJComp-355 379_28-A" => "CD3", 
    "FJComp-355 560_40-A" => "CD8", 

    "FJComp-355 820_60-A" => "CD4",
    "FJComp-355 670_30-A" => "CD4",

    "FJComp-640 780_60-A" => "CCR7",
    "FJComp-405 780_60-A" => "CD45RA", 

    "FJComp-561 780_60-A" => "CD127", 
    "FJComp-640 670_30-A" => "CD25", 

    "FJComp-561 610_20-A" => "Helios", 
    "FJComp-561 585_15-A" => "Foxp3", 
    "Foxp3-IgM" => "Foxp3",

    "FJComp-405 710_40-A" => "PD-1", 
    "FJComp-640 730_35-A" => "CXCR5", 

    "FJComp-405 670_30-A" => "CCR6", 
    "FJComp-488 715_30-A" => "CXCR3", 

    "FJComp-405 605_40-A" => "CCR4", 
    "FJComp-488 525_50-A" => "CCR10", 

    "FJComp-405 450_50-A" => "CD103", 
    "FJComp-355 740_35-A" => "CD69",
    "FJComp-405 515_20-A" => "HLA-DR"
])

data, labels, groups, gating = FlowWorkspace.load( files;
    workspace = workspace, channelMap = channelMap)

select!( labels, Not([ "CD4","CD8","Memory","Th17 | Th22","CD127- CD25+","non-Tregs"]))
select!( labels, Not(filter(name -> occursin("threshold", name), names(labels))))

###############################################################################
########################################################### calculate embedding
clusters, embedding = embed(data,path = "data/workspace.som",
    perplexity = 300, maxIter = 10000)

embeddingCoordinates = map(SVector, embedding[2,:], -embedding[1,:])
(xmin, xmax), (ymin, ymax) = extrema(embedding, dims = 2)
codes = select(hcat(labels, groups), AsTable(:) => ByRow(encode ∘ values) => "encoding")[:,:encoding]

###############################################################################
######################################################### colormaps and filters

##################################################### initialise color settings
channelRange, nlevels = range(-3, 7, length = 50), 10
toIndex(x::AbstractVector{<:Number}) = toIndex(x, channelRange; nlevels = nlevels)

colorIndex = combine(data, [ col => toIndex => col for col ∈ names(data) ])
channelLevels = range(extrema(channelRange)...,length=nlevels)

labelPalette = OrderedDict([ name=>'#'*hex(get(ColorSchemes.flag_gu,x)) for (name,x) ∈ zip(names(labels),range(0,1,length=size(labels,2))) ])
labelColors = channelview(fill(parse(RGBA{Float64}, "#DDDDDD"), size(data, 1)))
for (name,color) ∈ labelPalette colors!( Dict([ "name"=>name, "color"=>color ]) ) end

channelPalette = channelview(map(x -> RGBA(get(reverse(ColorSchemes.matter), x)), range(0, 1, length = nlevels)))
channelHexcodes = map(x->'#'*hex(x),colorview(RGB,view(channelPalette,1:3,:)))

##################################################### initalise filter settings
selection = OrderedDict([ name => true for name ∈ [names(labels);names(groups)] ])
selections = fill(true, size(data, 1))

populationNames = names(labels)
conditionNames = filter(name -> ~occursin(r"\d", name), names(groups))
groupNames = filter(name -> occursin(r"\d", name), names(groups))

###############################################################################
############################################################## construct canvas

const ol = JSServe.Dependency( :ol, # OpenLayers
    map( path -> joinpath(@__DIR__, path), [

        "assets/ol/ol.js",
        "assets/ol/ol.css",
        "assets/ol/sidebar.css"
    ])
)

const d3 = JSServe.Dependency( :d3, [ # data-driven documents
        joinpath(@__DIR__, "assets/d3/d3.v6.min.js"),
])

const style = JSServe.Dependency( :style, [ # custom styling todo(@gszep) remove remote dep
        "//cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.2/css/all.min.css",
        joinpath(@__DIR__, "assets/style.css"),
])

############################################## extensions loaded at end of body
extensions = map( extension -> HTML("""
    <script src="$(JSServe.Asset(joinpath(@__DIR__, extension)))" ></script>
"""), ["assets/ol/sidebar.js","assets/ol/colorbar.js","assets/Sortable.js","assets/violins.js","assets/boxplots.js","assets/utils.js"])

########################################################################## body
app = App() do session::Session
    Map = DOM.div(id = "map", class = "sidebar-map")
    
    JSServe.onload(session, Map, olMap([(xmin, xmax), (ymin, ymax)]; port = port))
    return DOM.div(id = "application", style, DOM.title("FlowAtlas.jl"), sidebar(session), Map, extensions)
end

###############################################################################
############################################# open app as locally hosted server
try
    global server = JSServe.Server(app, "127.0.0.1", port;
        verbose = true, routes = Routes( "/" => app,
            r"/assetserver/" * r"[\da-f]"^40 * r"-.*" => file_server,

            r"/\d+/\d+/\d+.png" => context -> tile(context;extrema = [(xmin, xmax),(ymin, ymax)]),
            r"/colors" => colors!, r"/selection" => selection!, r"/gate" => gate, r"/count" => count,

            r"/favicon.ico" => context -> HTTP.Response(500),
            r".*" => context -> response_404() )
    )

    browser(url)
    printstyled("[FlowAtlas] ", color = :blue)

    println("Running at $url ~ Happy exploring! ")
    wait()
        
catch exception
    if exception isa InterruptException
        close(server)

        printstyled("[FlowAtlas] ", color = :blue)
        println("Closing... have a lovely day! 🎈")

    else 
        rethrow()
    end
end