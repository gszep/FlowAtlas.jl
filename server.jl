import Pkg
Pkg.activate(".")

printstyled(color=:cyan,"[pkg] ")
printstyled("Installing/updating julia package dependencies...\n")

Pkg.instantiate()
Pkg.update()

printstyled("[FlowAtlas] ", color = :blue)
println("Compling Julia code... please be patient :) Will release pre-built sysimages in future")

using JSServe: Session, @js_str, file_server, response_404, Routes, App, DOM, Observable, on
using JSServe, HTTP, JSON, Images, FileIO, ImageIO
using GigaScatter, ColorSchemes

using Serialization: serialize,deserialize
using FlowWorkspace: inpolygon

using FlowWorkspace, StaticArrays, DataFrames, Glob
using StatsBase, GigaSOM, TSne
using StatsBase: normalize 

include("lib/components.jl")
include("lib/embed.jl")
include("lib/markers.jl")

include("lib/sys.jl")
include("lib/tile.jl")
include("lib/gate.jl")
include("lib/map.jl")

#################### server defaults
port = 3141
url = "http://localhost:$port"

#################### data inputs
workspace = "data/workspace.wsp"
files = glob"data/390C/*_BM_*.cleaned.fcs"
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

embeddingCoordinates = map( SVector, embedding[2,:], -embedding[1,:] )
(xmin,xmax), (ymin,ymax) = extrema(embedding,dims=2)

###############################################################################
######################################################## interactive components
fluorescence, automatic = Observable(true), Observable(true)

################################################## colormaps
channelRange, nlevels = range(-3, 7, length = 50), 10
toIndex(x::AbstractVector{<:Number}) = toIndex(x, channelRange;nlevels = nlevels)
colorIndex = combine(data, [ col => toIndex => col for col ∈ names(data) ])

segments = channelview(fill(parse(RGBA, "#EEEEEE00"), size(data, 1)))
palette = channelview( map( x->RGBA(get( reverse(ColorSchemes.curl),x)), range(0,1,length=nlevels)))

channel = Observable(first(names(data)))
colors = Observable(palette[:,colorIndex[!,channel[]] ])

legend = map( x-> RGB(get( ColorSchemes.Accent_8,x)), range(0,1,length=size(labels,2)))

legend = [ map( name -> Group(
    name,"#" * hex(popfirst!(legend))), names(labels) );
    map(name -> Group(name, true), names(groups))
]

############################# filter codes
codes = select( hcat(labels, groups), 
    AsTable(:) => ByRow(encode ∘ values) => "encoding")[:,:encoding]

labelCode = encode([i ≤ size(labels, 2) for i ∈ 1:length(legend)])
groupCode = encode([i > size(labels, 2) for i ∈ 1:length(legend)])

############################# update on channel change
on(channel) do channel
    if fluorescence[]
        
        selected = codes .& encode(map(x -> x.selected[], legend))
        selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
        
        markerColors = palette[ :, colorIndex[!,channel] ]
        markerColors[ :, .~selected ] .= 0.0
        colors[] = markerColors
    end
end

on(fluorescence) do fluorescence

    selected = codes .& encode(map(x -> x.selected[], legend))
    selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
    
    markerColors = fluorescence ? palette[ :, colorIndex[!,channel[]] ] : segments
    markerColors[ :, .~selected] .= 0.0
    colors[] = markerColors
end

############################# update legend interactions
for label ∈ legend
    
    mask = label.name ∈ names(labels) ? labels[!,label.name] : groups[!,label.name]
    if label.name ∈ names(labels) segments[:,mask] .= channelview([parse(RGBA, label.color[])]) end

    on(label.color) do labelColor
        segments[:,mask] .= channelview([parse(RGBA, labelColor)])
        
        selected = codes .& encode(map(x -> x.selected[], legend))
        selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
        
        markerColors = fluorescence[] ? palette[ :, colorIndex[!,channel[]] ] : segments
        markerColors[ :, .~selected] .= 0.0
        colors[] = markerColors
    end
    
    on(label.selected) do selected
        if selected segments[:,mask] .= channelview([parse(RGBA, label.color[])]) end
        
        selected = codes .& encode(map(x -> x.selected[], legend))
        selected = @. ( selected & labelCode ≠ 0 ) & ( selected & groupCode ≠ 0 )
        
        markerColors = fluorescence[] ? palette[ :, colorIndex[!,channel[]] ] : segments
        markerColors[ :, .~selected] .= 0.0
        colors[] = markerColors
    end
end

###############################################################################
############################################################## construct canvas

########################################################################## head
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
"""), ["assets/ol/sidebar.js","assets/Sortable.js","assets/violins.js","assets/uuid.js"])

########################################################################## body
app = App() do session::Session
    Map = DOM.div(id = "map", class = "sidebar-map")

    JSServe.onload( session, Map, olMap( [(xmin,xmax), (ymin,ymax)]; port=port) )
    return DOM.div( id="application", style, DOM.title("FlowAtlas.jl"), sidebar(session), Map, extensions )
end

###############################################################################
############################################# open app as locally hosted server
try
    global server = JSServe.Server(app, "127.0.0.1", port;
        verbose = true, routes = Routes( "/" => app,
            r"/assetserver/" * r"[\da-f]"^40 * r"-.*" => file_server,

            r"/\d+/\d+/\d+.png" => context -> tile(context;extrema=[(xmin,xmax),(ymin,ymax)]),
            r"/gate" => context -> gate(context),

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