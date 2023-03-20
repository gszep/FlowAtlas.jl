function colorbar(colorScale, height = 100, width = 100, origin = { x: 0, y: 0 },
    margin = { top: 5, right: 60, bottom: 5, left: 0 }) {

    function chart(containers) {
        containers.each(function (data) {
            var colorbar = d3.select(this).append("svg").classed("colorbar", true)

                .attr("x", origin.x - margin.right)
                .attr("y", origin.y - margin.top)

                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)

                .append("g").classed("colorbar", true)
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")")

            colorbar.append("g").attr("class", "axis")
            colorbar.append("g").attr("class", "bar")

            var bar = d3.range(0, height, by = height / (colorScale.domain().length - 1))
                .concat(height).reverse()

            var barScale = colorScale.copy()
            barScale.range(bar)

            /////////////////////////////////////////////// numberic axis
            colorbar.select(".axis")
                .attr("transform", "translate (" + width + "," + 0 + ")")
                .call(d3.axisRight(barScale))

            /////////////////////////////////////////////// fill in colormap
            colorbar.select(".bar").selectAll("colormap")
                .data(d3.range(0, height)).enter().append("rect")
                .attr("width", width).attr("height", 2)

                .style("fill", y => colorScale(barScale.invert(y)))
                .attr("y", y => y)
        })
    }

    return chart
}