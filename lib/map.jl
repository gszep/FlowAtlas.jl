function olMap(extrema::Array{<:Tuple}; port::Int = 3141)
    (xmin, xmax), (ymin, ymax) = extrema

    return js"""
        function (container){

            var colorScale = $d3.scaleLinear().domain($channelLevels).range($channelHexcodes)
            $d3.select("#map").call( colorbar( colorScale, height=300, width=20, origin={x:40,y:10} ))

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
                        return 'http://localhost:$port/'+[z,y,x].join("/")+'.png?'+'channel='+document.getElementById('channel').value+'&seed='+Math.random()
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

            document.querySelectorAll('#polygon-interactions button').forEach(
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

                                    boxplots( $labelCounts, new Set($populationNames), new Set($conditionNames),
                                        barcolors = { "CD4 | EM": "#afe3bd" },
                                        markercolors = { "390C": "#a6cee3", "403C": "#e31a1c", "412C": "#6a3d9a" })

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
                                selectedFeature ? overlay.setPosition(selectedFeature.getGeometry().getLastCoordinate()): overlay.setPosition(undefined)
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

    ############################################################# sortable legend
    Legend = DOM.div( id = "legend", class = "legend", map( group ->
                DOM.ul( id = group, DOM.h2( class = "legend-header", selected=true, uppercasefirst(group),

                    ############################################# de/select all in group
                    onclick = js"""
                        var element = document.getElementById($group)
                        var selected = element.querySelector("h2").getAttribute('selected') == 'true' ? 'false' : 'true'
                        Array.from(element.children).forEach( x => x.setAttribute('selected',selected) )

                        fetch( 'http://localhost:$port/selection', {

                            headers: { 'Content-Type': 'application/json' },
                            method: 'POST',

                            body: JSON.stringify({
                                selected: selected == 'true' ? true : false,
                                name: $group
                            })

                        }).then( response => {
                            document.getElementById('map').tiles.refresh()

                        }).catch( error => {
                            console.error('Error:',error)
                        })

                    """
                )),
            ["populations","conditions","groups"]
        )
    )

    JSServe.onload( session, Legend, js"""
        function (container){

            ////////////////////////////////////////////// populate legend with groups and event listeners
            var colors = $labelPalette
            for ( const [key,value] of Object.entries($([names(labels);names(groups)])) ) {

                var group = document.createElement('li')
                group.setAttribute('id',value)

                /////////////////////////////////////////////////// selection event
                group.setAttribute('selected',true)
                group.onclick = function (event) {
                    if (event.target.tagName != 'INPUT') {

                        var element = document.getElementById(value)
                        element.getAttribute('selected') == 'true' ? element.setAttribute('selected','false') : element.setAttribute('selected','true')
                        fetch( 'http://localhost:$port/selection', {

                            headers: { 'Content-Type': 'application/json' },
                            method: 'POST',

                            body: JSON.stringify({
                                selected: element.getAttribute('selected') == 'true' ? true : false,
                                name: value
                            })

                        }).then( response => {
                            document.getElementById('map').tiles.refresh()

                        }).catch( error => {
                            console.error('Error:',error)
                        })
                    }
                }

                var label = document.createElement('label')
                label.innerHTML = value

                var color = document.createElement('input')
                Object.assign( color, {

                    value: key < $(size(labels,2)) ? colors[key] : '#DDDDDD',
                    type: 'color',

                    ////////////////////////////////////////////// marker color event
                    onchange: function (event) {
                        fetch( 'http://localhost:$port/colors', {

                            headers: { 'Content-Type': 'application/json' },
                            method: 'POST',

                            body: JSON.stringify({
                                color: event.target.value,
                                name: value
                            })

                        }).then( response => {
                            document.getElementById('map').tiles.refresh()

                        }).catch( error => {
                            console.error('Error:',error)
                        })
                    }
                })

                group.appendChild(color)
                group.appendChild(label)

                if ( key < $(size(labels,2)) )
                    document.getElementById('populations').appendChild(group)

                else if ( value.match(/\d+/g) != null)
                    document.getElementById('groups').appendChild(group)

                else
                    document.getElementById('conditions').appendChild(group)
            }

            //////////////////////////////////////////////////// legend groups as sortable lists
            var legend = Object.assign(...['populations','conditions','groups'].map(
                group => ({

                    [group]: Sortable.create(
                        document.getElementById(group), {
                        group: group == 'populations' ? 'populations' : 'shared',

                        dataIdAttr: 'id',
                        animation: 150,

                        /////////////////////////////////////////// re-allocation of groups event
                        onAdd: function (event) {

                            // update_obs($populationNames,legend.populations)
                            // update_obs($conditionNames, legend.conditions)
                            // update_obs($groupNames,legend.groups)

                            fetch( 'http://localhost:$port/selection', {

                                headers: { 'Content-Type': 'application/json' },
                                method: 'POST',
                                body: JSON.stringify({ selected: NaN, name: '' })

                            }).then( response => {
                                document.getElementById('map').tiles.refresh()

                            }).catch( error => {
                                console.error('Error:',error)
                            })
                        }
                    })
                })
            ))
        }
    """)

    ############################################################# violins
    Violins = DOM.div(id = "violins", class = "svg-container")
    JSServe.onload( session, Violins, js"""
        function (container){

            var svg = d3.select("div#violins").append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")

                .attr("viewBox", "0 0 300 $(30*length(names(data)))")
                .classed("svg-content", true)
        }
    """)

    ############################################################# boxplots
    Boxplots = DOM.div(id = "boxplots", class = "svg-container")
    JSServe.onload( session, Boxplots, js"""
        function (container){

            var svg = d3.select("div#boxplots").append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")

                .attr("viewBox", "0 0 300 $(300*length(populationNames))")
                .classed("svg-content", true)
        }
    """)

    return DOM.div( id = "sidebar", class = "sidebar collapsed",

        ######################################## sidebar layout
        DOM.div(class = "sidebar-tabs",
            DOM.ul(role = "tablist", map( (href, class) ->
                HTML("""<li><a href="#$href" role="tab"><i class="$class"></i></a></li>"""),
                ["annotations","expression","frequency"], ["fab fa-amilia","fa fa-chart-area","fa fa-chart-pie"]
            )),

            DOM.ul(role = "tablist", map( (href, class) ->
                HTML("""<li><a href="#$href" role="tab"><i class="$class"></i></a></li>"""),
                ["settings"], ["fa fa-cog"]
            )),
        ),

        ######################################## sidebar widgets
        DOM.div(class = "sidebar-content",

            DOM.div(class = "sidebar-pane",id = "annotations",
                HTML("""<h1 class="sidebar-header">Annotations<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                DOM.div( id = "channel-selector", class="option-group",
        
                    "Colour by",
                    DOM.select( DOM.option.(["Labels";names(data)]), id="channel", onchange = js"""
                        document.querySelector('svg.colorbar').style.visibility = document.getElementById('channel').value == 'Labels' ? 'hidden' : 'visible'
                        document.getElementById("map").tiles.refresh()
                    """)
                ),

                Legend
            ),

            DOM.div(class = "sidebar-pane",id = "expression",
                HTML("""<h1 class="sidebar-header">Expression<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                html"""
                <div id="js-overlay" style="display:block">
                    <button id="js-remove">Remove</button>
                </div>
                <div id="polygon-interactions" class="option-group">
                    <button id="pan" type="button" class="primary">Pan</button>
                    <button id="polygon" type="button" class="success">Polygon</button>
                    <button id="modify" type="button" class="primary">Modify</button>
                    <button id="delete" type="button" class="danger">Delete</button>
                </div>
                """,
                Violins
            ),

            DOM.div(class = "sidebar-pane",id = "frequency",
                HTML("""<h1 class="sidebar-header">Frequency<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                Boxplots
            ),

            DOM.div(class = "sidebar-pane",id = "settings",
                HTML("""<h1 class="sidebar-header">Settings<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                Boxplots
            )
        )
    )
end