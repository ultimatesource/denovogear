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

    var activeNode = null;

    var format = d3.format(",.6e");

    d3.select("svg").remove();

    var zoom = d3.zoom()
      .on("zoom", zoomed);

    var chartWrapper = d3.select("#chart_wrapper");

    var svg = chartWrapper.append("svg");
    var width = $("svg").parent().width();
    var height = $("svg").parent().height();

    svg
        .attr("width", width)
        .attr("height", height);

    var container = svg.call(zoom)
      .append("g");

    var treesHidden = true;
    d3.select("#sample_tree_toggle").on("click", function() {
      var buttonSelection = d3.select(this);

      d3.selectAll(".sampleTree")
          .attr("visibility", function(d) {
            var ret;

            if (treesHidden) {
              buttonSelection
                  .attr("class", "btn btn-danger")
                  .text("Hide Trees");
              return "visible";
            }
            else {
              buttonSelection
                  .attr("class", "btn btn-success")
                  .text("Show Trees");
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
      .attr("dx", 15)
      .attr("dy", 15)
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

        var dngData = null;

        if (d.dataNode.data.dngOutputData !== undefined) {
          dngData = d.dataNode.data.dngOutputData;
        }
        else {
          // TODO: should probably be some sort of search for the correct
          // child, rather than assuming it's the first one.
          console.log(d.dataNode);
          dngData = d.dataNode.data.sampleIds.children[0].dngOutputData;
        }

        d3.select("#id_display").attr("value", d.dataNode.id);
        d3.select("#gt_display").attr("value", dngData.GT);
        d3.select("#gq_display").attr("value", dngData.GQ);
        d3.select("#gp_display").attr("value", dngData.GP);
        d3.select("#dp_display").attr("value", dngData.DP);
        d3.select("#mup_display").attr("value", format(dngData.MUP));
        d3.select("#mu1p_display").attr("value", format(dngData.MU1P));
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
