// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global utils */
/* global sampleTreeView */
/* exported PedigreeView */

var PedigreeView = (function(d3, PubSub) {
  "use strict";

  var format = d3.format(",.6e");

  var PedigreeView = function(parentElement, graphData) {
    this._graphData = graphData;
    this._parentElement = parentElement;
    this._create();
    this.update();

    PubSub.subscribe("DNG_OVERLAY_UPDATE", this._stateUpdate.bind(this));
    PubSub.subscribe("ACTIVE_NODE_CHANGED", this._stateUpdate.bind(this));
  };

  PedigreeView.prototype._create = function() {

    this._parentElement.append("svg")
        .attr("class", "pedigree-svg");

    var svg = this._parentElement.select(".pedigree-svg");

    svg.style("height", "100%").style("width", "100%");

    var zoom = d3.zoom()
      .on("zoom", zoomed.bind(this));
    this._container = svg.call(zoom)
      .append("g")
        .attr("class", "container");

    this._links_container = this._container.append("g")
        .attr("class", "links-container");

    this._nodes_container = this._container.append("g")
        .attr("class", "nodes-container");

    d3.select("#sample_tree_toggle").on("click", function() {
      PubSub.publish("SAMPLE_TREE_TOGGLE");
    });

    function zoomed() {
      this._container.attr("transform", d3.event.transform);
    }
  };
 

  PedigreeView.prototype.update = function() {

    this._updateLinks();
    this._updateNodes();

    PubSub.publish("PEDIGREE_VIEW_UPDATE", { activeNode: this._activeNode });
  };

  PedigreeView.prototype._updateLinks = function() {

    var links = this._graphData.links;

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
      .attr("dx", function(d) {
        return utils.halfwayBetween(d.source.x, d.target.x) - 45;
      })
      .attr("dy", function(d) {
        return utils.halfwayBetween(d.source.y, d.target.y);
      });

    visualLinksEnterUpdate.select("text")
      .text(function(d) {
        if (linkHasData(d)) {
          return d.dataLink.data.mutation;
        }
      });

    visualLinksEnterUpdate.select("path")
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

    function linkHasData(d) {
      return d.dataLink !== undefined && d.dataLink.data !== undefined;
    }
  };

  PedigreeView.prototype._updateNodes = function() {
    var visualNodesUpdate = this._nodes_container.selectAll(".node")
        .data(this._graphData.nodes);

    var visualNodesEnter = visualNodesUpdate.enter()
      .append("g")
        .attr("class", "node");

    var tree = sampleTreeView.createSampleTree();

    visualNodesEnter.selectAll(".yolo")
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

    visualNodesEnter.append("path")
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

    visualNodesEnter.append("text")
      .attr("dx", 15)
      .attr("dy", 15)
      .text(function(d) { 
        if (d.type !== "marriage") {
          return d.dataNode.id;
        }
      })
      .style("pointer-events", "none")
      .style("font", "10px sans-serif");

    visualNodesEnter.attr("transform", function(d) {
      return svgTranslateString(d.x, d.y);
    });
  };

  function nodeClicked(d) {
    if (d.type !== "marriage") {
      PubSub.publish("ACTIVE_NODE_CHANGED",
        { activeNode: d, activeNodeSelection: this });
    }
  }

  PedigreeView.prototype._stateUpdate = function(topic, data) {

    switch(topic) {

    case "DNG_OVERLAY_UPDATE":
      this.update();
      return;
    case "ACTIVE_NODE_CHANGED":
      this._activeNode = data.activeNode;
      d3.selectAll(".nodeSymbol").style("fill", fillColor);
      d3.select(data.activeNodeSelection).style("fill", "DarkSeaGreen");
      return;
    case "SAMPLE_TREE_TOGGLE":

      d3.selectAll(".sampleTree").attr("visibility", function() {
        if (state.showSampleTrees) {
          return "visible";
        }
        else {
          return "hidden";
        }
      });

      d3.select("#sample_tree_toggle")
        .attr("class", function() {
          if (state.showSampleTrees) {
            return "btn btn-danger";
          }
          else {
            return "btn btn-success";
          }
        })
        .text(function() {
          if (state.showSampleTrees) {
            return "Hide Trees";
          }
          else {
            return "Show Trees";
          }
        });
      return;
    default:
      console.log("unkown event");
      return;
    }
  };

  function svgTranslateString(x, y) {
    return "translate(" + x + "," + y + ")";
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

  return PedigreeView;

}(d3, PubSub));
