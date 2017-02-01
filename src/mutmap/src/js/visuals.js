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

    var treesHidden = true;

    d3.select("#sample_tree_toggle").on("click", function() {
      d3.selectAll(".sampleTree")
          .attr("visibility", function(d) {
            var ret;

            if (treesHidden) {
              return "visible";
            }
            else {
              return "hidden";
            }
          });

      treesHidden = !treesHidden;
    });

    var visualLinks = container
      .append("g")
        .attr("class", "links")
      .selectAll("line")
        .data(links)
      .enter().append("g")
        .attr("class", "link");

    visualLinks.append("path")
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
        .attr("stroke-width", function(d) {
          if (linkHasData(d)) {
            return 5;
          }
          return 1;
        })
        .attr("stroke", function(d) {
          if (linkHasData(d)) {
            return "green";
          }

          return "#999";
        })
        .attr("stroke-dasharray", function(d) {
          if (d.type == "duplicate") {
            return "5, 5";
          }
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
          return d.dataLink.data.mutation;
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


    var tree = createSampleTree();

    visualNodes.selectAll(".yolo")
        .data(function(d) {
          if (d.type == "person") {
            return [d];
          }
          else {
            return [null];
          }
        })
      .enter().each(function(d) {
        if (d) {
          var selection = d3.select(this);
          tree.node(d)(selection);
        }
      });

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
      .attr("visibility", function(d) {
        if (d.type === "marriage") {
          return "hidden";
        }
        else {
          return "visible";
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
      return svgTranslateString(d.x, d.y);
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

  function svgTranslateString(x, y) {
    return "translate(" + x + "," + y + ")";
  }

  return { doVisuals: doVisuals };

  function createSampleTree() {

    var node;
    var sampleTree;

    function my(selection) {

        var root = d3.hierarchy(node.dataNode.data.sampleIds);
        var treeWidth = 160;
        var treeHeight = 50;
        var cluster = d3.cluster().size([treeWidth, treeHeight]);
        cluster(root);

        var rootNode = root.descendants()[0];

        sampleTree = selection.append("g")
            .attr("class", "sampleTree")
            .attr("transform", function(d) {
              // center the tree with the tree's root node overlapping the
              // current node
              return svgTranslateString(-rootNode.x, -rootNode.y);
            })
            .attr("visibility", "hidden");

        var sampleTreeLinks = sampleTree.selectAll("sampleTreeLink")
            .data(root.descendants().slice(1))
          .enter().append("g")
            .attr("class", "sampleTreeLink")
          .append("line")
            .attr("stroke", "#999")
            .attr("x1", function(d) { return d.x; })
            .attr("y1", function(d) { return d.y; })
            .attr("x2", function(d) { return d.parent.x; })
            .attr("y2", function(d) { return d.parent.y; });

        var sampleTreeNodes = sampleTree.selectAll("sampleTreeNode")
            .data(root.descendants().slice(1))
          .enter().append("g")
            .attr("class", "sampleTreeNode")
            .attr("transform", function(d) {
              return svgTranslateString(d.x, d.y);
            });
          
        sampleTreeNodes.append("circle")
            .attr("r", 5)
            .attr("fill", "green");

        sampleTreeNodes.append("text")
            .attr("dx", 10)
            .attr("dy", ".35em")
            .text(function(d) {
              return d.data.name;
            });
    }

    my.node = function(value) {
      if (!arguments.length) return node;
      node = value;
      return my;
    };

    my.sampleTreeSelection = function() {
      return sampleTree;
    };

    return my;
  }

}());
