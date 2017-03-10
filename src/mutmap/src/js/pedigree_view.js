// eslint exceptions
//
/* global d3 */
/* global utils */
/* exported PedigreeView */

var PedigreeView = (function($, d3, store) {
  "use strict";

  store.subscribe(stateChanged);

  d3.select("#sample_tree_toggle").on("click", function() {
    var action = {
      type: "TOGGLE_SAMPLE_TREES"
    };
    store.dispatch(action);
  });

  var PedigreeView = function(selection, graphData) {

    this.nodes = graphData.nodes;
    this.links = graphData.links;

    var zoom = d3.zoom()
      .on("zoom", this.zoomed.bind(this));

    var svg = selection.append("svg").attr("class", "pedigree-svg");
    var width = $(".pedigree-svg").parent().width();
    var height = $(".pedigree-svg").parent().height();

    svg.attr("width", width).attr("height", height);

    this.container = svg.call(zoom)
      .append("g")
        .attr("class", "pedigree-container");

    this.linkContainer = this.container
      .append("g")
        .attr("class", "link-container");

    this.nodeContainer = this.container
      .append("g")
        .attr("class", "node-container");

    this.update();
  };


  PedigreeView.prototype.zoomed = function() {
      this.container.attr("transform", d3.event.transform);
  }

  PedigreeView.prototype.update = function() {

    var visualLinks = this.linkContainer 
      .selectAll(".link")
        .data(this.links);

    // update selection
    visualLinks.selectAll(".link__line")
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
        });

    // update selection
    visualLinks.selectAll(".link__text")
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

    var linksEnter = visualLinks.enter().append("g")
        .attr("class", "link");

    linksEnter.append("path")
        .attr("class", "link__line")
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
        .attr("y2", function(d) { return d.target.y; })
        .attr("stroke-width", function(d) { return 1; })
        .attr("stroke", function(d) { return "#999"; });

    linksEnter.append("text")
        .attr("class", "link__text")
        .attr("dx", function(d) {
          return utils.halfwayBetween(d.source.x, d.target.x) - 45;
        })
        .attr("dy", function(d) {
          return utils.halfwayBetween(d.source.y, d.target.y);
        });

    function linkHasData(d) {
      return d.dataLink !== undefined && d.dataLink.data !== undefined;
    }

    var visualNodes = this.nodeContainer.selectAll(".node")
      .data(this.nodes)
      .enter()
      .append("g")
        .attr("class", "node");

    var tree = sampleTreeView.createSampleTree();

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

    var nodeSymbols = visualNodes.append("path")
      .attr("class", "nodeSymbol")
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
      .style("fill", fillColor)
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
      return utils.svgTranslateString(d.x, d.y);
    });

  }

  function nodeClicked(d) {
    if (d.type !== "marriage") {
      var action = {
        "type": "NODE_CLICKED",
        "selection": this,
        "node": d
      };
      store.dispatch(action);
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

  function stateChanged() {

    var state = store.getState();

    if (state.activeNode) {
      d3.selectAll(".nodeSymbol").style("fill", fillColor);
      d3.select(state.activeNodeSelection).style("fill", "DarkSeaGreen");
    }

    d3.selectAll(".sampleTree").attr("visibility", function(d) {
      if (state.showSampleTrees) {
        return "visible";
      }
      else {
        return "hidden";
      }
    });

    d3.select("#sample_tree_toggle")
      .attr("class", function(d) {
        if (state.showSampleTrees) {
          return "btn btn-danger";
        }
        else {
          return "btn btn-success";
        }
      })
      .text(function(d) {
        if (state.showSampleTrees) {
          return "Hide Trees";
        }
        else {
          return "Show Trees";
        }
      });;
  }


  return PedigreeView;

}($, d3, store));
