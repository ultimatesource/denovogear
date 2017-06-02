var mutationDistributionView = (function(d3, PubSub, utils) {
  "use strict";

  var optionsManager = utils.createOptionsManager();

  function MutationDistributionView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'graphData', 'vcfText'],
      providedOptions: options
    });

    var distProc = new DistributionProcessor(options.vcfText);


    this._graphData = options.graphData;

    this._selection = options.renderInto;

    var container = this._selection.append("div")
        .attr("class", "row")
      .append("div")
        .attr("class", "mutation-dist-container col-xs-12");

    var svg = container.append("svg")
        .style("width", "100%")
        .style("height", "100%");

    var g = svg.append("g");

    var boundingRect = this._selection.node().getBoundingClientRect();
    var width = boundingRect.width;
    var height = boundingRect.height;

    this._links_container = g.append("g")
        .attr("class", "links-container");

    this._nodes_container = g.append("g")
        .attr("class", "nodes-container");

    this._updateLinks();
    this._updateNodes();
  }

  MutationDistributionView.prototype._updateLinks = function() {

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

          return "#999";
        });

    function linkHasMutation(d) {
      return d.dataLink !== undefined && d.dataLink.data !== undefined &&
        d.dataLink.data.mutation !== undefined;
    }

  MutationDistributionView.prototype._updateNodes = function() {

    var visualNodesUpdate = this._nodes_container.selectAll(".node")
        .data(this._graphData.nodes);

    var visualNodesEnter = visualNodesUpdate.enter()
      .append("g")
        .attr("class", "node");

    var visualNodesEnterUpdate = visualNodesEnter.merge(visualNodesUpdate);

    var tree = sampleTreeView.createSampleTree();

    visualNodesEnterUpdate.each(function(d) {

      var symbolSize = 500;

      if (d.type === "person") {
        if (d.dataNode.sex === "male") {
          d3.select(this).append("path")
              .attr("d", d3.symbol().type(d3.symbolSquare).size(symbolSize))
              .attr("fill", "steelblue");
        }
        else {
          d3.select(this).append("path")
              .attr("d", d3.symbol().type(d3.symbolCircle).size(symbolSize))
              .attr("fill", "tomato");
        }
      }
    });

    //visualNodesEnterUpdate.call(gpNode());
    //
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
      return utils.svgTranslateString(d.x, d.y);
    });
  };



  }
  
  function createMutationDistributionView(options) {
    return new MutationDistributionView(options);
  }


  function DistributionProcessor(vcfText) {

    //var vcfText;

    //fetch("platinum_filtered.vcf")
    //.then(function(response) {
    //  return response.text();
    //})
    //.then(function(text) {
    //  vcfText = text;

    //  processData();
    //});
    
    processData();

    function processData() {
      var startTime = new Date();

      var count = 0;
      var counts = {};
      //vcfParser.parseVCFText(vcfText);
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

          if (index % 1000 === 0) {
            console.log(dnl);
          }
        }
      });

      console.log(count);
      console.log(counts);

      var endTime = new Date();
      var elapsed = endTime - startTime;

      console.log("Parsing time:", elapsed / 1000);
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


  }

  return {
    createMutationDistributionView: createMutationDistributionView
  };

}(d3, PubSub, utils));
