function olMap(extrema::Array{<:Tuple}; port::Int=3141)
    (xmin,xmax), (ymin,ymax) = extrema

    return js"""
        function (container){

            ////////////////////////////////////////////////// tile layer
            var projection = new $ol.proj.Projection({
                code: 'raster',

                extent: $([ymin,-xmax,ymax,-xmin]),
                units: 'pixels',
            })

            $ol.proj.addProjection(projection);
            var tiles = new $ol.layer.Tile({

                source: new $ol.source.XYZ({
                    tileUrlFunction: function (zxy) {

                        [z,x,y] = zxy
                        return 'http://localhost:$port/'+[z,y,x].join("/")+'.png&seed='+Math.random()
                    },

                    projection: 'raster',
                    wrapX: false,
                })
            })

            ////////////////////////////////////////////////// polygon layer
            var polygons = new $ol.source.Vector({wrapX:false})
            var gates = new $ol.layer.Vector({
                source: polygons,

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

            //////////////////////////////////////////////////////// initialise map
            var map = new $ol.Map({
                layers: [ tiles, gates ],
                target: 'map',

                view: new $ol.View({
                    projection: 'raster',
                    center: [0, 0],
                    zoom: 0,
                }),
            })

            var overlay = new $ol.Overlay({ element: document.getElementById('js-overlay') })
            map.addOverlay(overlay) // delete button pop-up

            document.getElementById("map").tiles = tiles.getSource()
            document.getElementById("map").gates = gates.getSource()

            ///////////////////////////////////////////////////////// polygon interactions
            var selectedFeature
            var interaction
            var features

            document.querySelectorAll('#button-interactions button').forEach(
                button => button.addEventListener('click', function(event) {
                    var features

                    map.removeInteraction(interaction)
                    switch (event.target.id) {

                        case 'polygon':
                            interaction = new $ol.interaction.Draw({
                                type: 'Polygon', source: polygons })

                            map.addInteraction(interaction)
                            polygons.on('addfeature', function(event) {

                                console.assert( event.feature.getGeometry().getType()=="Polygon" )
                                event.feature.setId( uuidv4() )

                                ////////////////////////////////////////////// send gate to julia
                                fetch( 'http://localhost:$port/gate', {
                                    headers: { 'Content-Type': 'application/json' },
                                    method: 'POST',

                                    body: JSON.stringify({
                                        id: event.feature.getId(),
                                        coordinates: event.feature.getGeometry().getCoordinates()[0]
                                    })

                                //////////////////////////////////////////// update violins/stats
                                }).then( response => response.text() ).then( data => {
                                    var data = $d3.csvParse(data,$d3.autoType)

                                    violins({

                                        binCentres: Array.from( $d3.group(data,d=>d.channelValue).keys() ),
                                        density: Object.keys(data[0]).filter(x=>x!="channelValue").map( function(channel) {
      
                                            return {name: channel, values: data.map(x=>x[channel]).map(
                                                y=>y/$d3.max(data.map(x=>x[channel])) ) }
                                        })
                                    })

                                }).catch( error => {
                                    console.error('Error:',error)
                                })
                            })
                            break

                        case 'modify':
                            interaction = new $ol.interaction.Modify({
                                features: new $ol.Collection(polygons.getFeatures()) })
                            
                            map.addInteraction(interaction)
                            break

                        case 'delete':
                            interaction = new $ol.interaction.Select({
                                condition: $ol.events.condition.click, layers: [gates] })

                            map.addInteraction(interaction)
                            interaction.on('select', function(event) {

                                selectedFeature = event.selected[0]
                                selectedFeature ? overlay.setPosition(selectedFeature.getGeometry().getExtent()): overlay.setPosition(undefined)
                            })
                            break

                        default:
                            break
                    }

                    var snap = new $ol.interaction.Snap({source: polygons })
                    map.addInteraction(snap)
                })
            )

            document.getElementById('js-remove').addEventListener('click', function() {
                polygons.removeFeature(selectedFeature)

                overlay.setPosition(undefined)
                interaction.getFeatures().clear()
            })

            var sidebar = new $ol.contr$ol.Sidebar({ element: 'sidebar', position: 'left' })
            map.addControl(sidebar)
        }
    """
end

function sidebar(session::Session)

    ############################################################# interactive legend
    Legend = DOM.div(id="legend",class="legend", HTML("""

        <ul id="populations">
            <li>foo 1</li>
            <li>foo 2</li>
            <li>foo 3</li>
        </ul>

        <ul id="conditions">
            <li>foo 1</li>
            <li>foo 2</li>
            <li>foo 3</li>
        </ul>

        <ul id="groups">
            <li>bar 1</li>
            <li>bar 2</li>
            <li>bar 3</li>
        </ul>    
    """)
    )

    JSServe.onload( session, Legend, js"""
        function (container){

            var populations = document.getElementById('populations')
            Sortable.create(populations, {

                group: 'populations',
                animation: 100
            })

            var conditions = document.getElementById('conditions')
            Sortable.create(conditions, {

                group: 'shared',
                animation: 100
            })

            var groups = document.getElementById('groups')
            Sortable.create(groups, {

                group: 'shared',
                animation: 100
            })
        }
    """)

    ############################################################# violins
    Violins = DOM.div(id = "violins", class="svg-container")
    JSServe.onload( session, Violins, js"""
        function (container){

            var svg = d3.select("div#violins").append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")

                .attr("viewBox", "0 0 300 $(30*length(names(data)))")
                .classed("svg-content", true)
        }
    """)

    return DOM.div( id = "sidebar", class = "sidebar collapsed",

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
                DOM.div( class="container",
                        
                    DOM.label("Population"),
                    DOM.label(class="switch",
                            
                        DOM.input(type="checkbox",checked=true,
                        onchange=js"""update_obs($fluorescence,this.checked);document.getElementById("map").tiles.refresh()"""),
                        DOM.span(class="slider round")
                    ),
        
                    DOM.label("Fluorescence Intensity"),
                    DOM.select(DOM.option.(names(data)),
                    onchange=js"""update_obs($channel,this.options[this.selectedIndex].text);document.getElementById("map").tiles.refresh()""")
                ),

                Legend,
                DOM.div( class = "container",
    
                    ######################################### interactive legend
                    DOM.span( DOM.label(class = "switch", style = "width:70px",	
                        DOM.input(type = "checkbox", checked = true,
                            onchange = js"""$(map(legend->legend.selected,
                                filter( label->label.name ∈ names(labels), legend))
                                ).map(x=>{update_obs(x,this.checked);document.getElementById("map").tiles.refresh()})"""
                        ),
                                
                        DOM.span(class = "slider text", style = "font-weight:bold", "Populations")),
                        map( legend -> DOM.div( class = "container",
                                    
                            DOM.input(type = "color",name = legend.name,id = legend.name,
                                value = legend.color, onchange = js"""
                                    update_obs($(legend.color),this.value);document.getElementById("map").tiles.refresh()"""),
                            DOM.label(class = "switch", style = "width:100%",
                                        
                            DOM.input(type = "checkbox",
                                checked = legend.selected, onchange = js"""
                                    update_obs($(legend.selected),this.checked);document.getElementById("map").tiles.refresh()"""),
                            DOM.span(class = "slider text", legend.name))
        
                        ), filter(label -> label.name ∈ names(labels), legend) ),
                    ),
                        
                    DOM.span( DOM.label(class = "switch", style = "width:70px",	
                        DOM.input(type = "checkbox", checked = true, 					
                            onchange = js"""$(map(legend->legend.selected,
                                filter( label->label.name ∈ names(groups), legend))
                                ).map(x=>{update_obs(x,this.checked);document.getElementById("map").tiles.refresh()})"""
                        ),
                                
                        DOM.span(class = "slider text", style = "font-weight:bold", "Groups")),
                        map( legend -> DOM.div( class = "container",
                                    
                            DOM.input(type = "color",name = legend.name,id = legend.name,
                                value = legend.color, onchange = js"""
                                    update_obs($(legend.color),this.value);document.getElementById("map").tiles.refresh()"""),
                            DOM.label(class = "switch", style = "width:100%",
                                        
                            DOM.input(type = "checkbox",
                                checked = legend.selected, onchange = js"""
                                    update_obs($(legend.selected),this.checked);document.getElementById("map").tiles.refresh()"""),
                            DOM.span(class = "slider text", legend.name))
        
                        ), filter(label -> label.name ∈ names(groups), legend) ),
                    ),
                )
            ),

            DOM.div(class = "sidebar-pane",id = "comparisons",
                HTML("""<h1 class="sidebar-header">Comparisons<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                html"""
                <div id="js-overlay" style="display:block">
                    <button id="js-remove">Remove</button>
                </div>
                <div id="button-interactions" class="btn-group btn-group-sm" role="group" aria-label="Draw">
                    <button id="pan" type="button" class="btn btn-primary">Pan</button>
                    <button id="polygon" type="button" class="btn btn-success">Polygon</button>
                    <button id="modify" type="button" class="btn btn-primary">Modify</button>
                    <button id="delete" type="button" class="btn btn-danger">Delete</button>
                </div>
                """,
                Violins
            ),

            html"""
            <div class="sidebar-pane" id="clustering">
                <h1 class="sidebar-header">Clustering<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>
            </div>

            <div class="sidebar-pane" id="settings">
                <h1 class="sidebar-header">Settings<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>
            </div>"""
        )
    )
end