// eslint exceptions
//
/* global d3 */
/* global utils */
/* global sampleTreeView */

var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  var format = d3.format(",.6e");

  mutmap.PedigreeView = Backbone.View.extend({
    
    initialize: function() {

      this.d3el = d3.select(this.el);

      this.activeNodeModel = this.model.get('activeNodeModel');

      this.model.on('change', this.render.bind(this));
      this.activeNodeModel.on('change', this._updateActiveNode.bind(this));

      this._showSampleTrees = false;

      this.d3el.append("svg")
          .attr("class", "pedigree-svg");

      this.dim = utils.getDimensions(this.d3el);
      this._svg = this.d3el.select(".pedigree-svg")
          .attr("width", this.dim.width)
          .attr("height", this.dim.height);

      this._container = this._svg.append("g");

      this._links_container = this._container.append("g")
          .attr("class", "links-container");

      this._nodes_container = this._container.append("g")
          .attr("class", "nodes-container");

      this._gpNode = gpNode();


      this.render();

      // Center and zoom the pedigree
      // bounding box of the links in the graph serves as a heuristic for
      // calculating the center of the graph
      var boundingBox = this._links_container.node().getBBox();
      var centerX = boundingBox.width / 2;
      var centerY = boundingBox.height / 2;

      var xCorrection = (this.dim.width / 2) - centerX;
      var yCorrection = (this.dim.height / 2) - centerY;

      var zoom = d3.zoom()
        .on("zoom", zoomed.bind(this));
      this._svg.call(zoom);

      zoom.translateBy(this._svg, xCorrection, yCorrection);
      zoom.scaleBy(this._svg, 0.1, 0.1);

      var transition = this._svg.transition().duration(750);
      zoom.scaleTo(transition, 1.0);
    },

    render: function() {

      this.dim = utils.getDimensions(this.d3el);
      this._svg
          .attr("width", this.dim.width)
          .attr("height", this.dim.height);

      this._updateLinks();
      this._updateNodes();
      this._updateSampleTrees();
      this._updateActiveNode();

    },

    _updateLinks: function() {

      var self = this;

      var links = this.model.get('graphData').links;

      var visualLinksUpdate = this._links_container.selectAll(".link")
          .data(links);

      var visualLinksEnter = visualLinksUpdate.enter()
        .append("g")
          .attr("class", "link");

      var visualLinksEnterUpdate = visualLinksEnter.merge(visualLinksUpdate);

      visualLinksEnter.append("path")
          .attr("d", function(d) {

            var points = self._genLinkPoints(d);

            if (d.type == "child" || d.type == "spouse") {
              return "M" + points.source.x + "," + points.source.y +
                "L" + points.target.x + "," + points.target.y;
            }
            else {

              var controlX = utils.halfwayBetween(points.source.x,
                points.target.x);
              // TODO: parameterize the control point Y value by the distance
              // between the nodes, rather than hard coding
              var controlY = points.source.y - 100;
              return "M" + points.source.x + "," + points.source.y
                + "Q" + controlX + "," + controlY + ","
                + points.target.x + "," + points.target.y;
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
          var points = self._genLinkPoints(d);
          return utils.halfwayBetween(points.source.x, points.target.x);
        })
        .attr("dy", function(d) {
          var points = self._genLinkPoints(d);
          return utils.halfwayBetween(points.source.y, points.target.y);
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
              return d3.schemeCategory20[7];
            }

            return "#999";
          });

      function linkHasMutation(d) {
        return d.dataLink !== undefined && d.dataLink.data !== undefined &&
          d.dataLink.data.mutation !== undefined;
      }
    },

    _updateNodes: function() {

      var self = this;
      var visualNodesUpdate = this._nodes_container.selectAll(".node")
          .data(this.model.get('graphData').nodes);

      var visualNodesEnter = visualNodesUpdate.enter()
        .append("g")
          .attr("class", "node");

      var visualNodesEnterUpdate = visualNodesEnter.merge(visualNodesUpdate);

      var tree = sampleTreeView.createSampleTree();

      // TODO: this is a hack. I don't even really know why it works...
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


      var nodeSymbols = visualNodesEnter.append("g")
          .attr("class", function(d) {
            if (d.type === "person") {
              return "node-symbol node-symbol__" + d.dataNode.sex;
            }
          })
          .on("click", this._nodeClicked.bind(this));

      var nodeSymbolsEnterUpdate = nodeSymbols.merge(visualNodesUpdate);

      nodeSymbolsEnterUpdate.call(this._gpNode);
      
      visualNodesEnter.append("text")
        .attr("x", 12)
        .attr("y", 18)
        .text(function(d) { 
          if (d.type !== "marriage") {
            return d.dataNode.id;
          }
        })
        .style("pointer-events", "none")
        .style("font", "8px sans-serif");

      visualNodesEnter.attr("transform", function(d) {

        var points = self._genNodePoints(d);
        return utils.svgTranslateString(points.x, points.y);
      });
    },

    _updateSampleTrees: function() {

      var showSampleTrees = this.model.get('showSampleTrees');

      d3.selectAll(".sampleTree")
          .attr("visibility", function() {
            if (showSampleTrees) {
              return "visible";
            }
            else {
              return "hidden";
            }
          });

      d3.select("#sample_tree_toggle")
          .attr("class", function() {
            if (showSampleTrees) {
              return "btn btn-danger";
            }
            else {
              return "btn btn-success";
            }
          })
          .text(function() {
            if (showSampleTrees) {
              return "Hide Trees";
            }
            else {
              return "Show Trees";
            }
          });


    },

    _updateActiveNode: function() {

      // Create variable to preserve 'this' context
      var activeNodeModel = this.activeNodeModel;
      d3.selectAll(".node-symbol").each(function(d) {
        if (d.dataNode === activeNodeModel.get('dataNode')) {
          d3.select(this).classed("node-symbol--selected", true);
        }
        else {
          d3.select(this).classed("node-symbol--selected", false);
        }
      });
    },

    _nodeClicked: function(d) {
      if (d.type !== "marriage") {
        this.activeNodeModel.set({
          dataNode: d.dataNode
        });
      }
    },

    // Converts x and y (which are ratios between 0 and 1) into relative
    // values based on the view dimensions
    _genNodePoints: function(d) {
      return this._genPoints(d.x, d.y);
    },

    _genLinkPoints: function(d) {
      return {
        source: this._genPoints(d.source.x, d.source.y),
        target: this._genPoints(d.target.x, d.target.y)
      }
    },

    _genPoints(xRatio, yRatio) {
      // Apply 20% padding
      var width = this.dim.width - (.2 * this.dim.width);
      return {
        x: xRatio * width,
        y: yRatio * (width / 3)
      };
    },

  });

  function zoomed() {
    this._container.attr("transform", d3.event.transform);
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

  function gpNode() {

    var nodeWidth = 20;

    var color = d3.scaleOrdinal(d3.schemeCategory10);
    var pie = d3.pie()
      .sort(null);
    var piePath = d3.arc()
      .outerRadius(nodeWidth/2)
      .innerRadius(0);

    function my(selection) {

      selection.each(function(d) {

        if (d.type !== "person") {
          return;
        }
        
        if (d.dataNode.sex === "male") {
          if (d.dataNode.data.dngOutputData) {

            var squareWidth = nodeWidth;
            var squareHeight = squareWidth;

            var partitions = calculatePartitions(
              d.dataNode.data.dngOutputData.GP, squareWidth);

            var squareUpdate = d3.select(this).selectAll('.partition')
                .data(partitions);

            var squareEnter = squareUpdate.enter()
              .append("rect")
                .attr("class", "partition")

            var squareEnterUpdate = squareEnter.merge(squareUpdate);

            squareEnterUpdate
                .attr("x", function(d) {
                  return -(squareWidth / 2) + d.offset;
                })
                .attr("y", -(squareHeight / 2))
                .attr("width", function(d) {
                  return d.split*squareWidth;
                })
                .attr("height", squareHeight)
                .attr("fill", function(d, i) {
                  return d3.schemeCategory10[i];
                });
          }
        }

        else if (d.dataNode.sex === "female") {
          if (d.dataNode.data.dngOutputData) {
            var gpSplits =
              calculateGpSplits(d.dataNode.data.dngOutputData.GP);

            var pieUpdate = d3.select(this).selectAll(".arc")
                .data(pie(gpSplits));

            var pieEnter = pieUpdate.enter().append("path")
                .attr("class", "arc");

            var pieEnterUpdate = pieEnter.merge(pieUpdate);

            pieEnterUpdate
                .attr("d", piePath)
                .attr("fill", function(d, i) { 
                  return d3.schemeCategory10[i];
                });
          }
        }
      });
    }

    return my;
  }

  function calculatePartitions(gp, length) {
    var gpSplits = calculateGpSplits(gp);
    var offsets = calculateOffsets(gpSplits, length);

    return gpSplits.map(function(split, index) {
      return {
        split: split,
        offset: offsets[index]
      };
    });
  }

  function calculateOffsets(splits, length) {

    var offsets = [];
    var offset = 0;
    splits.forEach(function(split, index) {
      offsets.push(offset);
      offset += split*length;
    });

    return offsets;
  }

  function calculateGpSplits(gp) {

    var gpIndicesTable;

    // TODO: Generate this table algorithmically.
    if (gp.length === 3) {
      // According to VCF spec, ordering is AA, AB, BB
      // see GL section in
      // https://samtools.github.io/hts-specs/VCFv4.2.pdf
      gpIndicesTable = [
        [ 1 ],
        [ 2 ]
      ];
    }
    else if (gp.length === 6) {
      // According to VCF spec, ordering is AA, AB, BB, AC, BC, CC
      // see GL section in
      // https://samtools.github.io/hts-specs/VCFv4.2.pdf
      gpIndicesTable = [
        [ 1, 3 ],
        [ 2, 4, 5 ]
      ];
    }
    else if (gp.length === 10) {
      gpIndicesTable = [
        [ 1, 3, 6 ],
        [ 2, 4, 5, 7, 8, 9 ]
      ];
    }
    else {
      throw "GP indices table not defined for GP length: " + gp.length;
    }

    var refOnly = gp[0];

    var includesRef = 0;
    gpIndicesTable[0].forEach(function(index) {
      includesRef += gp[index];
    });

    var doesNotIncludeRef = 0;
    gpIndicesTable[1].forEach(function(index) {
      doesNotIncludeRef += gp[index];
    });

    return [
      refOnly,
      includesRef,
      doesNotIncludeRef
    ];
  }


  function createPedigreeView(options) {
    return new PedigreeView(options);
  }

  return {
    createPedigreeView: createPedigreeView
  };

}(Backbone, d3, utils));
