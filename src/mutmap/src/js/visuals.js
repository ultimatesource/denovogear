// eslint exceptions
//
/* global d3 */
/* global utils */
/* exported visuals */

var visuals = (function() {
  "use strict";
  
  function doVisuals(graphData) {

    var nodes = graphData.nodes;
    var links = graphData.links;

    var height = 500;
    var activeNode = null;

    d3.select("svg").remove();

    var zoom = d3.zoom()
      .on("zoom", zoomed);


    var chartWrapper = d3.select("#chart_wrapper");
    var dim = chartWrapper.node().getBoundingClientRect();

    var svg = chartWrapper.append("svg")
        .attr("width", dim.width)
        .attr("height", height);

    var container = svg.call(zoom)
      .append("g");

      
    var visualLinks = container
      .append("g")
        .attr("class", "links")
      .selectAll("line")
      .data(links)
      .enter()
      .append("g")
        .attr("class", "link");

    visualLinks.append("line")
      .attr("stroke-width", function(d) {
        if (linkHasData(d)) {
          //return 5;
        }
        return 1;
      })
      .attr("stroke", function(d) {
        if (linkHasData(d)) {
          //return "green";
        }

        return "#999";
      })
      .attr("x1", function(d) { return d.source.x; })
      .attr("y1", function(d) { return d.source.y; })
      .attr("x2", function(d) { return d.target.x; })
      .attr("y2", function(d) { return d.target.y; });

    visualLinks.append("text")
      .attr("dx", function(d) {
        return utils.halfwayBetween(d.source.x, d.target.x) - 45;
      })
      .attr("dy", function(d) {
        return utils.halfwayBetween(d.source.y, d.target.y);
      })
      .text(function(d) {
        if (linkHasData(d)) {
          //return d.dataLink.data.mutation;
        }
      });

    function linkHasData(d) {
      return d.dataLink !== undefined && d.dataLink.data !== undefined;
    }


    var visualNodes = container
      .append("g")
        .attr("class", "nodes")
      .selectAll(".node")
      .data(nodes)
      .enter()
      .append("g")
        .attr("class", "node");

    visualNodes.append("path")
      .attr("d", d3.symbol()
        .type(function(d) {
          if (d.type === "person") {
            if (d.dataNode.sex === "male") {
              return d3.symbolSquare;
            }
            else if (d.dataNode.sex === "female") {
              return d3.symbolCircle;
            }
          }
          else {
            return d3.symbolTriangle;
          }
        })
        .size(500))
      .attr("fill", fillColor)
      .attr("opacity", function(d) {
        if (d.type === "marriage") {
          return 0;
        }
        else {
          return 1;
        }
      })
      .on("click", nodeClicked);

    visualNodes.append("text")
      .attr("dx", 20)
      .attr("dy", ".35em")
      .text(function(d) { 
        if (d.type !== "marriage") {
          return d.dataNode.id;
        }
      });

    visualNodes.attr("transform", function(d) {
      return "translate(" + d.x + "," + d.y + ")";
    });

    // Create sample ID tree diagram for each node
    nodes.forEach(function(node) {
      if (node.type == 'person') {
        var root = d3.hierarchy(node.dataNode.data.sampleIds);
        var cluster = d3.cluster().size([100, 100]);
        cluster(root);
        //console.log(node);
        node.tree = root;

        var className = "sampleTree" + String(node.dataNode.id);
        var nodeClass = "sampleTreeNode" + String(node.dataNode.id);
        console.log(root.descendants());
        var sampleTreeNodes = svg.append("g")
            .attr("class", className)
            .data(root.descendants())
          .enter().append("circle")
            .attr("cx", function(d) { return d.x; })
            .attr("cy", function(d) { return d.y; })
            .attr("r", 5)
            .attr("class", nodeClass);
        //console.log(className);
      }
    });

    function zoomed() {
      container.attr("transform", d3.event.transform);
    }

    function nodeClicked(d) {
      if (d.type !== "marriage") {
        d3.select(activeNode).style("fill", fillColor);
        activeNode = this;
        d3.select(this).style("fill", "DarkSeaGreen");

        if (d.dataNode.data.dngOutputData !== undefined) {
          document.getElementById("id_display").value =
            d.dataNode.id;
        }
        else {
          document.getElementById("id_display").value = "";
        }

      }
    }

    function fillColor(d) {
      if (d.type === "person") {
        if (d.dataNode.sex === "male") {
          return "SteelBlue";
        }
        else if (d.dataNode.sex === "female") {
          return "Tomato";
        }
      }
      else {
        return "black";
      }
    }
  }

  return { doVisuals: doVisuals };

}());
