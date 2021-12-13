function boxplots(data, populations, conditions, batches,
    barcolors = {}, markercolors = {},
    margin = ({ top: 10, right: 50, bottom: 30, left: 50 }), padding = 1/4 ) {

    const svg = d3.select("div#boxplots svg")
    svg.selectAll("*").remove()

    var color = d3.scaleOrdinal().domain( batches.size>0 ? batches : ["Ungrouped"]).range( batches.size>0 ? [...batches].map(x => markercolors[x]) : ['#663F46'] )
    var height = populations.size * 256
    var width = 256

    var populationName = d3.scaleBand().domain(populations)
        .range([margin.top, height])

    var conditionName = d3.scaleBand().domain(conditions)
        .range([margin.left, width]).padding(padding)

    const axes = svg.append("g").selectAll("g").data(populations).join("g")
        .attr("transform", name => `translate(0,${populationName.step()
            - margin.bottom + 1.1 * populationName(name)})`)

    axes.append("text").attr("transform", "rotate(-90)")
        .attr("y", 0).attr("x", populationName.step() / 2 + margin.top - margin.bottom)
        .attr("dy", "1em").style("text-anchor", "middle")

        .style("font-family", "sans-serif").style("font-size", "14px")
        .text(name => name)

    axes.append("g").call(
        d3.axisBottom(conditionName).tickSizeOuter(0)).call(

            g => g.selectAll("text").attr("y", 10).attr("x", 0).attr("dy", ".35em")
                .attr("transform", "rotate(-45)")).style("text-anchor", "end")

    axes.selectAll("legend")
        .data(batches).enter().call(g => {

            g.append("circle").attr("cx", width + 10)
                .attr("cy", function (d, i) { return i * 25 - 3 * populationName.step() / 4 - margin.top + margin.bottom })

                .attr("r", 2).style("fill", color)
                .attr("stroke", "black").attr("stroke-width", "0.25px")

            g.append("text").attr("x", width + 10)
                .attr("y", function (d, i) { return i * 25 - 3 * populationName.step() / 4 - margin.top + margin.bottom })
                .style("fill", color).text(x => x).attr("dx", "0.5em").attr("dy", "0.1em")

                .style("alignment-baseline", "middle")
                .style("font-family", "sans-serif").style("font-size", "12px")
        })

    var frequencies = [], quantiles = []
    populations.forEach(population => {
        conditions.forEach(condition => {

            if (batches.size > 0) {
                batches.forEach(group => {

                    var frequency = d3.sum(data,
                        x => [population, condition, group].every(y => x.names.includes(y)) ? x.count : 0

                    ) / d3.sum(data,
                        x => [condition, group].every(y => x.names.includes(y)) ? x.count : 0
                    )

                    if (!Number.isNaN(frequency))
                        frequencies.push({
                            frequency: frequency,

                            population: population,
                            condition: condition,
                            group: group

                        })
                })
            } else {
                var frequency = d3.sum(data,
                    x => [population, condition].every(y => x.names.includes(y)) ? x.count : 0

                ) / d3.sum(data,
                    x => [condition].every(y => x.names.includes(y)) ? x.count : 0
                )

                if (!Number.isNaN(frequency))
                    frequencies.push({
                        frequency: frequency,

                        population: population,
                        condition: condition,
                        group: "Ungrouped"

                    })
            }

            var points = frequencies.filter(x => x.population == population).filter(x => x.condition == condition)
            quantiles.push({

                upper: d3.quantile(points.map(x => x.frequency).sort(d3.ascending), 3 / 4)+1e-2*d3.max(points.map(x => x.frequency).sort(d3.ascending)),
                lower: d3.quantile(points.map(x => x.frequency).sort(d3.ascending), 1 / 4)-1e-2*d3.max(points.map(x => x.frequency).sort(d3.ascending)),

                max: d3.max(points.map(x => x.frequency).sort(d3.ascending)),
                min: d3.min(points.map(x => x.frequency).sort(d3.ascending)),

                population: population,
                condition: condition,
            })
        })
    })

    var frequency = {}
    axes.append("g").call(
        g => g.attr("transform", `translate(${margin.left},0)`).each(function (name) {

            var max = d3.max(frequencies.filter(z => z.population.includes(name)), y => y.frequency)
            var scale = d3.scaleLinear().domain([0, max]).range([0, margin.bottom - populationName.step()])

            frequency[name] = scale
            return d3.select(this).call(d3.axisLeft(scale).tickSizeOuter(0)
                .ticks(width / 80).tickFormat(x => `${100 * x.toFixed(2)}%`))
        })
    )

    svg.selectAll("quantiles").data(quantiles).enter().call(g => {

        g.append("rect")
            .attr("x", x => conditionName(x.condition))
            .attr("y", x => frequency[x.population](x.upper)
                + 1.1 * populationName(x.population) + populationName.step() - margin.bottom)

            .attr("height", x => frequency[x.population](x.lower) - frequency[x.population](x.upper))
            .attr("width", conditionName.bandwidth())

            .style("fill", x => barcolors[x.population])
            .style("opacity", 0.5).style("cursor","pointer")

            .on("click", function(event,data){
                console.log(data.population,data.condition)
            })

        g.append("line")
            .attr("stroke", x => barcolors[x.population])
            .attr("stroke-width", "3px")

            .attr("x1", x => conditionName(x.condition) + conditionName.bandwidth() / 2)
            .attr("x2", x => conditionName(x.condition) + conditionName.bandwidth() / 2)

            .attr("y1", x => frequency[x.population](x.max)
                + 1.1 * populationName(x.population) + populationName.step() - margin.bottom)
            .attr("y2", x => frequency[x.population](x.min)
                + 1.1 * populationName(x.population) + populationName.step() - margin.bottom)
    })

    svg.selectAll("points").data(frequencies).enter().append("circle")

        .attr("cx", x => conditionName(x.condition)
            + conditionName.bandwidth() * random(padding, 1 - padding))

        .attr("cy", x => frequency[x.population](x.frequency)
            + 1.1 * populationName(x.population) + populationName.step() - margin.bottom)

        .attr("r", 2).style("fill", x => color(x.group))

    return svg.node()
}