function sidebar(session::Session, names::NamedTuple, colors::NamedTuple; port::Int = 3141)

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
            var colors = $(collect(values(colors.labels.names)))
            for ( const [key,value] of Object.entries($([names.populations;names.conditions;names.groups])) ) {

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

                    value: key < $(length(names.populations)) ? colors[key] : '#663F46',
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

                if ( key < $(length(names.populations)) )
                    document.getElementById('populations').appendChild(group)

                else if ( value.match(/\d+/g) != null)
                    document.getElementById('groups').appendChild(group)

                else
                    document.getElementById('conditions').appendChild(group)
            }

            //////////////////////////////////////////////////// legend groups as sortable lists
            ['populations','conditions','groups'].map( group => {

                Sortable.create( document.getElementById(group), {
                    group: group == 'populations' ? 'populations' : 'shared',

                    dataIdAttr: 'id',
                    animation: 150,

                    /////////////////////////////////////////// re-allocation of groups event
                    onAdd: function (event) {
                        fetch( 'http://localhost:$port/selection', {

                            headers: { 'Content-Type': 'application/json' },
                            method: 'POST',

                            body: JSON.stringify({
                                name: {
                                    'populations': Array.from(document.getElementById('populations').querySelectorAll('li')).map(x=>x.id),
                                    'conditions': Array.from(document.getElementById('conditions').querySelectorAll('li')).map(x=>x.id),
                                    'groups': Array.from(document.getElementById('groups').querySelectorAll('li')).map(x=>x.id)
                                }
                            })

                        }).then( response => {
                            document.getElementById('map').tiles.refresh()

                        }).catch( error => {
                            console.error('Error:',error)
                        })
                    }
                })
            })
        }
    """)

    ############################################################# violins
    Violins = DOM.div(id = "violins", class = "svg-container")
    JSServe.onload( session, Violins, js"""
        function (container){

            var svg = $d3.select("div#violins").append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")

                .attr("viewBox", "0 0 300 $(30*length(names.channels))")
                .classed("svg-content", true)
        }
    """)

    ############################################################# boxplots
    Boxplots = DOM.div(id = "boxplots", class = "svg-container")
    JSServe.onload( session, Boxplots, js"""
        function (container){

            var svg = $d3.select("div#boxplots").append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")

                .attr("viewBox", "0 0 300 $(300*length(names.populations))")
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
                    DOM.select( DOM.option.(["Labels";names.channels]), id="channel", onchange = js"""
                        document.querySelector('svg.colorbar').style.visibility = document.getElementById('channel').value == 'Labels' ? 'hidden' : 'visible'
                        document.getElementById("map").tiles.refresh()
                    """)
                ),

                Legend
            ),

            DOM.div(class = "sidebar-pane",id = "expression",
                HTML("""<h1 class="sidebar-header">Expression<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),

                HTML(""" <div id="js-overlay" style="display:block"> <button id="js-remove">Remove</button></div>"""),
                DOM.div( id="polygon-interactions", class="option-group",

                    DOM.button( id="pan", type="button", class="primary", "Pan"),
                    DOM.button( id="polygon", type="button", class="success", "Polygon"),
                    DOM.button( id="modify", type="button", class="primary", style="visibility: hidden", "Modify"),
                    DOM.button( id="delete", type="button", class="danger", "Delete"),

                    DOM.button( id="save-boxplots", type="button", "Save",
                        onclick=js"""

                            var script = document.createElement('script');
                            Object.assign( script, { type: 'text/javascript', src: $FileSaver })
                            document.body.appendChild(script)

                            var svg = document.querySelector("#violins svg").innerHTML
                            var head = '<svg title="violins" version="1.1" xmlns="http://www.w3.org/2000/svg">'

                            var blob = new Blob([head+svg+"</svg>"],{type:"image/svg+xml"})
                            saveAs(blob,"violins.svg")
                    """)
                ),
                Violins
            ),

            DOM.div(class = "sidebar-pane",id = "frequency",
                HTML("""<h1 class="sidebar-header">Frequency<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),
                DOM.div( id = "frequency-options", class="option-group",

                    DOM.button( id="frequency-calculate", type="button", "Calculate",
                        onclick=js"""

                            fetch( 'http://localhost:$port/count', {
                                method: 'GET'

                            //////////////////////////////////////////// update stats
                            }).then( response => response.json() ).then( data => {

                                var populations = Array.from(document.getElementById('populations').querySelectorAll('li')).filter(x=>x.getAttribute('selected')=='true')
                                var conditions = Array.from(document.getElementById('conditions').querySelectorAll('li')).filter(x=>x.getAttribute('selected')=='true')
                                var groups = Array.from(document.getElementById('groups').querySelectorAll('li')).filter(x=>x.getAttribute('selected')=='true')

                                boxplots( data,

                                    new Set(populations.map(x=>x.id)),
                                    new Set(conditions.map(x=>x.id)),
                                    new Set(groups.map(x=>x.id)),

                                    barcolors = populations.length > 0 ? Object.assign( ...populations.map( x => ({[x.id]: x.querySelector('input').value}) )) : {},
                                    markercolors = groups.length > 0 ? Object.assign( ...groups.map( x => ({[x.id]: x.querySelector('input').value}) )) : {}
                                )

                            }).catch( error => {
                                console.error('Error:',error)
                            })
                        """
                    ),

                    DOM.button( id="save-boxplots", type="button", "Save",
                        onclick=js"""

                            var script = document.createElement('script');
                            Object.assign( script, { type: 'text/javascript', src: $FileSaver })
                            document.body.appendChild(script)

                            var svg = document.querySelector("#boxplots svg").innerHTML
                            var head = '<svg title="boxplots" version="1.1" xmlns="http://www.w3.org/2000/svg">'

                            var blob = new Blob([head+svg+"</svg>"],{type:"image/svg+xml"})
                            saveAs(blob,"boxplots.svg")
                    """)
                ),
                Boxplots
            ),

            DOM.div(class = "sidebar-pane",id = "settings",
                HTML("""<h1 class="sidebar-header">Settings<span class="sidebar-close"><i class="fa fa-caret-left"></i></span></h1>"""),

                HTML("""<h2 align="center">Marker Size</h2>"""),
                DOM.input( type="range", min="0.01", step="0.01", max="5", value="1", class="slider", id="scale", onchange = js"""
                    document.getElementById("map").tiles.refresh()
                """)
            )
        )
    )
end