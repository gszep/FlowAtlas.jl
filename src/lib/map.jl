function olMap(extrema::Array{<:Tuple}, colors::NamedTuple; port::Int = 3141)
    (xmin, xmax), (ymin, ymax) = extrema

    return js"""
        function (container){

            var colorScale = $d3.scaleLinear().domain($(colors.channels.levels)).range($(colors.channels.hex))
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
                    fill: new $ol.style.Fill({ color: 'rgba(0, 136, 170, 0.2)' }),
                    stroke: new $ol.style.Stroke({ color: 'rgb(0, 136, 170)', width: 3 }),
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
            document.getElementById("map").gates = gates

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
      
                                            return { id: event.feature.getId(),
                                                name: channel, values: data.map(x=>x[channel]).map(
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

                var id = CSS.escape(selectedFeature.getId())
                document.querySelectorAll("#"+id).forEach( violin=> violin.remove())

                overlay.setPosition(undefined)
                interaction.getFeatures().clear()
            })

            var sidebar = new $ol.control.Sidebar({ element: 'sidebar', position: 'left' })
            map.addControl(sidebar)
        }
    """
end