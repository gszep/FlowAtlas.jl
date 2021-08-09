"""
Start FlowAtlas using:
```julia
julia> FlowAtlas.run()
```
"""
module FlowAtlas
include_dependency("../Project.toml")

using JSServe: Session, @js_str, file_server, response_404, Routes, App, DOM, Observable, on
using JSServe, HTTP, JSON, Images, FileIO, ImageIO
using GigaScatter, ColorSchemes

using Serialization: serialize,deserialize
using FlowWorkspace: inpolygon
using Base: NamedTuple

using FlowWorkspace, StaticArrays, DataFrames, MetaGraphs, OrderedCollections
using GigaSOM, TSne, StatsBase, Impute

include("lib/colors.jl")
include("lib/selection.jl")
include("lib/embed.jl")

include("lib/tile.jl")
include("lib/gate.jl")
include("lib/count.jl")

include("lib/sys.jl")
include("lib/map.jl")
include("lib/sidebar.jl")

############################################## javascript dependencies
const ol = JSServe.Dependency( :ol, [# OpenLayers
    joinpath(@__DIR__,"assets/ol/ol.js"), joinpath(@__DIR__,"assets/ol/ol.css")
])

const d3 = JSServe.Dependency( :d3, [ # data-driven documents
    joinpath(@__DIR__, "assets/d3/d3.v6.min.js"),
])

const FileSaver = JSServe.Dependency( :FileSaver, [ # saving files
    joinpath(@__DIR__, "assets/FileSaver.min.js")
])

const style = JSServe.Dependency( :style, [ # custom styling todo(@gszep) remove remote dep
    "//cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.2/css/all.min.css", joinpath(@__DIR__, "assets/style.css"),
])

const extensions = JSServe.Dependency( :extensions, map(  extension -> joinpath(@__DIR__,extension),
    ["assets/ol/sidebar.css","assets/ol/colorbar.js","assets/Sortable.js","assets/violins.js","assets/boxplots.js","assets/utils.js"]
))

function run( path::String; files::String=joinpath(dirname(path),"*.fcs"), transform::Function=x->asinh(x/250),
        port::Int = 3141, url::String = "http://localhost:$port", cols::Symbol=:union, drop::Union{Vector{String},Nothing}=nothing,
        nlevels::Int=10, channelRange = range(-3,7,length=50), channelScheme=reverse(ColorSchemes.matter), labelScheme=ColorSchemes.seaborn_colorblind,
        perplexity=300, maxIter=10000 )

    indexTransform(x::AbstractVector{<:Union{Number,Missing}}) = toIndex(x, channelRange; nlevels=nlevels)

    @info "Loading FCS files..."
    data, labels, groups, gating = FlowWorkspace.load( path; files=files, transform=transform, cols=cols)

    keep = map( graph -> map( idx -> MetaGraphs.get_prop(graph,idx,:name), filter(idx -> MetaGraphs.outdegree(graph,idx) == 0, 1:MetaGraphs.nv(graph))), values(gating) )
    select!( labels, isnothing(drop) ? union(keep...,["Unlabelled"]) : Not(drop) )

    names = (
        channels = Base.names(data), populations = Base.names(labels),
        conditions = filter(name -> ~occursin(r"\d", name), Base.names(groups)),
        groups = filter(name -> occursin(r"\d", name), Base.names(groups))
    )

    ###############################################################################
    ########################################################### calculate embedding
    som_path, _ = splitext(path)
    _, embedding = embed(data, path = som_path*".som", perplexity = perplexity, maxIter = maxIter)
    embedding = (
        coordinates = map(SVector, embedding[2,:], -embedding[1,:]),
        array = embedding, extrema = extrema(embedding,dims=2)
    )

    ###############################################################################
    ######################################################### colormaps and filters
    colors = (
        labels = (
            names = OrderedDict([ name=>'#'*hex(get(labelScheme,x)) for (name,x) ∈ zip(Base.names(labels),range(0,1,length=size(labels,2))) ]),
            rows = channelview(fill(parse(RGBA{Float64},"#663F46"), size(data,1))),
        
        ),
        channels = (
            levels = range(extrema(channelRange)...,length=nlevels),
            rows = combine(data, [ col => indexTransform => col for col ∈ Base.names(data) ]),
            colors = channelview([map(x->RGBA(get(channelScheme, x)),range(0,1,length=nlevels)); RGBA(0,0,0,0)]),
            hex = map(x->'#'*hex(x),colorview(RGB,view(channelview(map(x->RGBA(get(channelScheme, x)),range(0,1,length=nlevels))),1:3,:)))
            
        )
    )

    for (name,color) ∈ colors.labels.names colors!( Dict(["name"=>name,"color"=>color]),labels,groups,colors) end

    ##################################################### initalise filter settings
    selections = (
        codes = select(hcat(labels,groups), AsTable(:) => ByRow(encode ∘ values) => "encoding")[:,:encoding],
        names = OrderedDict([ name => true for name ∈ [Base.names(labels);Base.names(groups)] ]),
        rows = fill(true, size(data, 1))
    )

    ########################################################################## body
    app = App() do session::Session
        Map = DOM.div(id = "map", class = "sidebar-map")
        
        JSServe.onload(session, Map, olMap( embedding.extrema, colors; port=port))
        return DOM.div(id = "application", extensions, style, DOM.title("FlowAtlas.jl"), Map,

            HTML("""<script src="$(JSServe.Asset(joinpath(@__DIR__,"assets/ol/sidebar.js")))" ></script>"""),
            sidebar(session, names, colors; port=port)
        )
    end

    ###############################################################################
    ############################################# open app as locally hosted server
    try
        global server = JSServe.Server(app, "127.0.0.1", port;
            verbose = true, routes = Routes( "/" => app,
                r"/assetserver/" * r"[\da-f]"^40 * r"-.*" => file_server,

                r"/\d+/\d+/\d+.png" => x->tile(x,embedding,selections,colors,names),
                r"/colors" => x->colors!(x,labels,groups,colors), r"/selection" => x->selection!(x,selections,names),

                r"/gate" => x->gate(x,data,embedding,selections;channelRange=channelRange),
                r"/count" => x->count(x,selections),

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
end

@info """\n
    Welcome to FlowAtlas 🎈
    Start FlowAtlas server using:
        julia> FlowAtlas.run("workspace.wsp"; files="*.fcs")
\n"""
end