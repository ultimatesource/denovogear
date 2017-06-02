var mutationDistributionView = (function(d3, PubSub, utils) {
  "use strict";

  var optionsManager = utils.createOptionsManager();

  function MutationDistributionView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'graphData'],
      providedOptions: options
    });

    var graphData = options.graphData;

    var links = graphData.links;
    var nodes = graphData.nodes;
    console.log(graphData);

    this._selection = options.renderInto;

    var container = this._selection.append("div")
        .attr("class", "row")
      .append("div")
        .attr("class", "col-xs-12");

    var svg = container.append("svg")
        .style("width", "100%")
        .style("height", "100%");

    var g = svg.append("g");

    var boundingRect = this._selection.node().getBoundingClientRect();
    var width = boundingRect.width;
    var height = boundingRect.height;

    this._links_container = g.append("g")
        .attr("class", "links-container");

    var visualLinksUpdate = this._links_container.selectAll(".link")
        .data(links);

    var visualLinksEnter = visualLinksUpdate.enter()
      .append("g")
        .attr("class", "link");

    var visualLinksEnterUpdate = visualLinksEnter.merge(visualLinksUpdate);

    visualLinksEnter.append("path")
        .attr("d", function(d) {
          if (d.type == "child" || d.type == "spouse") {
            return "M" + d.source.x + "," + d.source.y +
              "L" + d.target.x + "," + d.target.y;
          }
          else {
            var controlX = utils.halfwayBetween(d.source.x, d.target.x);
            // TODO: parameterize the control point Y value by the distance
            // between the nodes, rather than hard coding
            var controlY = d.source.y - 100;
            return "M" + d.source.x + "," + d.source.y
              + "Q" + controlX + "," + controlY + ","
              + d.target.x + "," + d.target.y;
          }
        })
        .attr("fill", "transparent")
        .attr("stroke-dasharray", function(d) {
          if (d.type == "duplicate") {
            return "5, 5";
          }
        })
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });


    visualLinksEnter.append("text")
      .attr("text-anchor", "middle")
      .attr("dx", function(d) {
        return utils.halfwayBetween(d.source.x, d.target.x);
      })
      .attr("dy", function(d) {
        return utils.halfwayBetween(d.source.y, d.target.y);
      });

    visualLinksEnterUpdate.select("text")
      .text(function(d) {
        if (linkHasMutation(d)) {
          return d.dataLink.data.mutation;
        }
      });

    visualLinksEnterUpdate.select("path")
        .attr("stroke-width", function(d) {
          if (linkHasMutation(d)) {
            return 5;
          }
          return 1;
        })
        .attr("stroke", function(d) {
          if (linkHasMutation(d)) {
            return "#5cc464";
          }

          return "#999";
        });

    function linkHasMutation(d) {
      return d.dataLink !== undefined && d.dataLink.data !== undefined &&
        d.dataLink.data.mutation !== undefined;
    }

  }
  
  function createMutationDistributionView(options) {
    return new MutationDistributionView(options);
  }

  return {
    createMutationDistributionView: createMutationDistributionView
  };

}(d3, PubSub, utils));
