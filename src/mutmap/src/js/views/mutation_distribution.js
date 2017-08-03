var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  mutmap.MutationDistributionView = Backbone.View.extend({

    initialize: function(options) {

      var distProc = new DistributionProcessor(options.vcfText);

      this._counts = distProc.getCounts();

      this._maxCount = 0;
      Object.keys(this._counts).forEach(function(key) {
        if (this._counts[key] > this._maxCount) {
          this._maxCount = this._counts[key];
        }
      }, this);
      this._barScale = d3.scaleLinear()
        .domain([0, this._maxCount])
        .range([0, 80]);

      this._graphData = options.graphData;

      this.d3el = d3.select(this.el);

      var container = this.d3el.append("div")
          .attr("class", "row")
        .append("div")
          .attr("class",
            "mutation-dist-container col-xs-12 panel panel-default");

      var svg = container.append("svg")
          .style("width", "100%")
          .style("height", "100%");

      var zoom = d3.zoom()
        .on("zoom", zoomed.bind(this));
      svg.call(zoom);

      var g = svg.append("g");

      var boundingRect = this.d3el.node().getBoundingClientRect();
      var width = boundingRect.width;
      var height = boundingRect.height;

      this._links_container = g.append("g")
          .attr("class", "links-container");

      this._nodes_container = g.append("g")
          .attr("class", "nodes-container");

      this.render();

      function zoomed() {
        g.attr("transform", d3.event.transform);
      }
    },

    render: function() {
      this._updateLinks();
      this._updateNodes();
    },

    _updateLinks: function() {

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

            return "#ccc";
          });

      function linkHasMutation(d) {
        return d.dataLink !== undefined && d.dataLink.data !== undefined &&
          d.dataLink.data.mutation !== undefined;
      }

    },

    _updateNodes: function() {

      var visualNodesUpdate = this._nodes_container.selectAll(".node")
          .data(this._graphData.nodes);

      var visualNodesEnter = visualNodesUpdate.enter()
        .append("g")
          .attr("class", "node")
          .attr("transform", function(d) {
            return utils.svgTranslateString(d.x, d.y);
          });

      //var visualNodesEnterUpdate = visualNodesEnter.merge(visualNodesUpdate);

      visualNodesEnter.each(function(d) {

        var symbolSize = 500;

        if (d.type === "person") {
          if (d.dataNode.sex === "male") {
            d3.select(this).append("path")
                .attr("d", d3.symbol().type(d3.symbolSquare).size(symbolSize))
                .attr("fill", d3.schemeCategory20[1]);
          }
          else {
            d3.select(this).append("path")
                .attr("d", d3.symbol().type(d3.symbolCircle).size(symbolSize))
                .attr("fill", d3.schemeCategory20[7]);
          }
        }
      });

      visualNodesEnter.append("text")
          .attr("x", 15)
          .attr("y", 15)
          .text(function(d) { 
            if (d.type !== "marriage") {
              return d.dataNode.id;
            }
          })
          .style("pointer-events", "none")
          .style("font", "10px sans-serif");

      var barWidth = 10;
      var counts = this._counts;
      var scale = this._barScale;
      visualNodesEnter.append("rect")
          .attr("x", -25)
          .attr("y", function(d) {
            if (counts[d.dataNode.id]) {
              return -scale(counts[d.dataNode.id]);
            }
            else {
              return 0;
            }
          })
          .attr("width", barWidth)
          .attr("height", function(d) {
            if (counts[d.dataNode.id]) {
              return scale(counts[d.dataNode.id]);
            }
            else {
              return 0;
            }
          })
          .attr("fill", d3.schemeCategory20[4]);

      visualNodesEnter.append("text")
          .attr("x", -25 + barWidth/2)
          .attr("y", function(d) {

            var offset;
            if (counts[d.dataNode.id]) {
              offset = -scale(counts[d.dataNode.id]);
            }
            else {
              offset = 0;
            }

            return offset - 3;
          })
          .attr("text-anchor", "middle")
          .text(function(d) {
            if (counts[d.dataNode.id]) {
              return counts[d.dataNode.id];
            }
            else {
              return 0;
            }
          });

    }

  });

  
  function createMutationDistributionView(options) {
    return new MutationDistributionView(options);
  }


  function DistributionProcessor(vcfText) {

    this.counts = processData();

    function processData() {
      //var startTime = new Date();

      var count = 0;
      var counts = {};
      var lines = vcfText.split("\n");
      lines.forEach(function(line, index) {
        if (!line.startsWith("#") && line.length > 0) {

          var data = parseDataLine(line);
          var info = parseInfoColumn(data[7]);
          var dnl = parseDNL(info);

          if (counts[dnl.value] === undefined) {
            counts[dnl.value] = 0;
          }
          counts[dnl.value]++;
          count++;

          //if (index % 1000 === 0) {
          //  console.log(dnl);
          //}
        }
      });

      //var endTime = new Date();
      //var elapsed = endTime - startTime;
      //console.log("Parsing time:", elapsed / 1000);

      var processedCounts = {};

      Object.keys(counts).forEach(function(key) {
        var id = idFromKey(key);

        if (id.length > 0) {
          if (processedCounts[id] === undefined) {
            processedCounts[id] = 0;
          }

          processedCounts[id] += counts[key];
        }
      });

      return processedCounts;
    }

    function parseDataLine(line) {
      var columns = line.split("\t");
      return columns;
    }

    function parseInfoColumn(column) {
      return column.split(';');
    }

    function parseDNL(info) {
      return parsePair(info[6]);
    }

    function parsePair(pair) {
      var pairArray = pair.split('=');
      return {
        key: pairArray[0],
        value: pairArray[1]
      };
    }

    function idFromKey(key) {
      var idStartIndex = key.indexOf("/") + 1;
      var idStopIndex = key.indexOf("_");
      if (idStopIndex === -1) {
        idStopIndex = key.length;
      }

      return key.slice(idStartIndex, idStopIndex);
    }


  }

  DistributionProcessor.prototype.getCounts = function() {
    return this.counts;
  };

  return {
    createMutationDistributionView: createMutationDistributionView
  };

}(Backbone, d3, utils));
