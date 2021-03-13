function violins(data, color="#0088AA55", margin=({top:40,right:20,bottom:40,left:60}) ) {
  const svg = d3.select("div#violins svg")

  var height = data.density.length * 30
  var width = 300
  
  var groupName = d3.scalePoint().domain( data.density.map(d=>d.name) )
    .range([margin.top, height-margin.bottom])
  
  var binCentre = d3.scaleLinear().domain( d3.extent(data.binCentres) )
    .range([margin.left,width-margin.right])
  
  var density = d3.scaleLinear().domain( [0,1] )
    .range([0,-groupName.step()/2])
  
  svg.append("g").call(
    g => g.attr("transform", `translate(0,${20+height-margin.bottom})`).call(
      d3.axisBottom(binCentre).ticks(width/100) ))
  
  svg.append("g").call(
    g => g.attr("transform", `translate(${margin.left},0)`).call(
      d3.axisLeft(groupName).tickSize(0).tickPadding(10)).call(
      g=>g.select(".domain").remove() ))
  
  const group = svg.append("g").selectAll("g").data(data.density)
    .join("g").attr("transform", y=>`translate(0,${groupName(y.name)})`)

  var violin = d3.area().curve(d3.curveBasis)
    .x( (x,i)=> binCentre(data.binCentres[i]) )
    .y0( y=>-density(y) ).y1( y=>density(y) )
  
  group.append("path").attr("fill",color).attr("stroke",color)
    .attr("d", x=>violin(x.values) )
  return svg.node()
}