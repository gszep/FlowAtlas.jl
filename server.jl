import Pkg
Pkg.activate(".")

printstyled(color=:cyan,"[pkg] ")
printstyled("Installing/updating julia package dependencies...\n")

Pkg.instantiate()
Pkg.update()

printstyled("[FlowAtlas] ", color = :blue)
println("Compling Julia code... please be patient :) Will release pre-built sysimages in future")

using JSServe: @js_str, file_server, response_404, Routes, App, DOM, Observable, on
using JSServe, HTTP, GigaScatter, Plots, Images, FileIO, ImageIO

using Serialization: serialize,deserialize
using FlowWorkspace, StaticArrays, DataFrames, Glob
using GigaSOM, TSne

include("lib/components.jl")
include("lib/embed.jl")
include("lib/markers.jl")

include("lib/sys.jl")
include("lib/tile.jl")

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

###############################################################################
######################################################## interactive components
doneCalculations = Observable(true)
fluorescence, automatic = Observable(true), Observable(true)
compareGates = JSServe.Button("Compare Gates")

################################################## colormaps
channelRange, nlevels = range(-2, 7, length = 50), 10
toIndex(x::AbstractVector{<:Number}) = toIndex(x, channelRange;nlevels = nlevels)
colorIndex = combine(data, [ col => toIndex => col for col ∈ names(data) ])

segments = channelview(fill(parse(RGBA, "#EEEEEE00"), size(data, 1)))
palette = channelview(cgrad(:curl, nlevels, categorical = true, rev = true).colors.colors)

channel = Observable(first(names(data)))
colors = Observable(palette[:,colorIndex[!,channel[]] ])

legend = convert( Vector{RGB},
    cgrad(:Accent_8, size(labels, 2), categorical = true).colors.colors)

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

const style = JSServe.Dependency( :style, [ # custom styling
        "//cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.2/css/all.min.css",
        joinpath(@__DIR__, "assets/style.css"),
])

############################################## extensions loaded at end of body
extensions = map( extension -> HTML("""
    <script src="$(JSServe.Asset(joinpath(@__DIR__, extension)))" ></script>
"""), ["assets/ol/sidebar.js"])

app = App() do session::Session

    ###################################################### mapview construction
    Map = DOM.div(id = "map", class = "sidebar-map")
    JSServe.onload( session, Map, js"""
        function (container){

            var tiles = new $ol.layer.Tile({
                source: new $ol.source.XYZ({

                    url: 'http://localhost:$port/{z}/{y}/{x}.png',
                    wrapX: false,
                })
            })

            var source = new $ol.source.Vector()
            var gates = new $ol.layer.Vector({
                source: source,

                style: new $ol.style.Style({
                    fill: new $ol.style.Fill({
                        color: 'rgba(255, 255, 255, 0.2)',
                    }),
                    stroke: new $ol.style.Stroke({
                        color: '#ffcc33',
                        width: 2,
                    }),
                    image: new $ol.style.Circle({
                        radius: 7,
                        fill: new $ol.style.Fill({
                            color: '#ffcc33',
                        }),
                    }),
                }),
            })

            var map = new $ol.Map({
                layers: [ tiles, gates ],
                target: 'map',

                view: new $ol.View({
                    center: [0, 0],
                    zoom: 0,
                }),
            })

            document.getElementById("map").tiles = tiles
            var sidebar = new $ol.control.Sidebar({ element: 'sidebar', position: 'left' })
            map.addControl(sidebar)

            ///////////////////////////////////////////////////////// gating interactions
            var modify = new $ol.interaction.Modify({source: source});
            map.addInteraction(modify);

            var draw, snap; // global so we can remove them later
            var typeSelect = 'Polygon';

            function addInteractions() {
                draw = new $ol.interaction.Draw({
                    source: source,
                    type: 'Polygon',
                })
                map.addInteraction(draw)
                snap = new $ol.interaction.Snap({source: source})
                map.addInteraction(snap)
            }

            /**
             * Handle change event.
             */
            // typeSelect.onchange = function () {
            //     map.removeInteraction(draw)
            //     map.removeInteraction(snap)
            //     addInteractions()
            // }

            addInteractions()
        }
    """)

    return DOM.div( id = "application", style,
        DOM.div( id = "sidebar", class = "sidebar collapsed",

            ######################################## sidebar layout
            DOM.div(class = "sidebar-tabs",
                DOM.ul(role = "tablist", map( (href, class) ->
                    HTML("""<li><a href="#$href" role="tab"><i class="$class"></i></a></li>"""),
                    ["annotations","comparisons","clustering"], ["fab fa-amilia","fa fa-adjust","fa fa-arrows-alt"]
                )),

                DOM.ul(role = "tablist", map( (href, class) ->
                    HTML("""<li><a href="#$href" role="tab"><i class="$class"></i></a></li>"""),
                    ["settings"], ["fa fa-gear"]
                )),
            ),

            ######################################## sidebar widgets
            DOM.div(class = "sidebar-content",

                DOM.div(class = "sidebar-pane",id = "annotations",
                    HTML("""<h1 class="sidebar-header">Annotations<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                    DOM.div( class = "container",
        
                        ######################################### interactive legend
                        DOM.span( DOM.label(class = "switch", style = "width:70px",	
                            DOM.input(type = "checkbox", checked = true,
                                onchange = js"""$(map(legend->legend.selected,
                                    filter( label->label.name ∈ names(labels), legend))
                                    ).map(x=>JSServe.update_obs(x,this.checked))"""
                            ),
                                    
                            DOM.span(class = "slider text", style = "font-weight:bold", "Populations")),
                            map( legend -> DOM.div( class = "container",
                                        
                                DOM.input(type = "color",name = legend.name,id = legend.name,
                                    value = legend.color, onchange = js"""
                                        JSServe.update_obs($(legend.color),this.value)"""),
                                DOM.label(class = "switch", style = "width:100%",
                                            
                                DOM.input(type = "checkbox",
                                    checked = legend.selected, onchange = js"""
                                        JSServe.update_obs($(legend.selected),this.checked)"""),
                                DOM.span(class = "slider text", legend.name))
            
                            ), filter(label -> label.name ∈ names(labels), legend) ),
                        ),
                            
                        DOM.span( DOM.label(class = "switch", style = "width:70px",	
                            DOM.input(type = "checkbox", checked = true, 					
                                onchange = js"""$(map(legend->legend.selected,
                                    filter( label->label.name ∈ names(groups), legend))
                                    ).map(x=>JSServe.update_obs(x,this.checked))"""
                            ),
                                    
                            DOM.span(class = "slider text", style = "font-weight:bold", "Groups")),
                            map( legend -> DOM.div( class = "container",
                                        
                                DOM.input(type = "color",name = legend.name,id = legend.name,
                                    value = legend.color, onchange = js"""
                                        JSServe.update_obs($(legend.color),this.value)"""),
                                DOM.label(class = "switch", style = "width:100%",
                                            
                                DOM.input(type = "checkbox",
                                    checked = legend.selected, onchange = js"""
                                        JSServe.update_obs($(legend.selected),this.checked)"""),
                                DOM.span(class = "slider text", legend.name))
            
                            ), filter(label -> label.name ∈ names(groups), legend) ),
                        ),
                    )
                ),

                html"""
                <div class="sidebar-pane" id="comparisons">
                    <h1 class="sidebar-header">Comparisons<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>
                </div>

                <div class="sidebar-pane" id="clustering">
                    <h1 class="sidebar-header">Clustering<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>
                </div>

                <div class="sidebar-pane" id="settings">
                    <h1 class="sidebar-header">Settings<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>
                </div>"""
            )
        ),
    Map, extensions )
end

###############################################################################
    ############################################# open app as locally hosted server
try
    global server = JSServe.Server(app, "127.0.0.1", port;
        verbose = true, routes = Routes( "/" => app,
            r"/assetserver/" * r"[\da-f]"^40 * r"-.*" => file_server,

            r"/\d+/\d+/\d+.png" => context -> tile(context;extrema = extrema(embedding, dims = 2)),
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