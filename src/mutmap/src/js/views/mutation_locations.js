
var mutationLocationsView = (function(d3, PubSub, utils) {
  "use strict";

  function MutationLocationsView(options) {

    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'mutationLocationData'],
      providedOptions: options
    });

    var parent = options.renderInto;

    var container = parent.append("div")
        .attr("class", "row")
      .append("div")
        .attr("class",
              "mutation-location-container col-xs-12 panel panel-default")

    var svg = container.append("svg")
        .style("width", "100%")
        .style("height", "100%");

    var boundingRect = parent.node().getBoundingClientRect();
    var width = boundingRect.width;
    var height = boundingRect.height;

    var margins = {
      left: 20,
      top: 20,
      right: 20,
      bottom: 20
    };

    var chromWidth = width - (margins.left + margins.right);

    var g = svg.append("g")
        .attr("transform",
              utils.svgTranslateString(margins.left, margins.top));
    
    //g.call(sampleMutationList().data(options.mutationLocationData));

    SampleMutationListView.create({
      renderInto: g,
      data: options.mutationLocationData[0],
      width: chromWidth,
    });

  }

  function SampleMutationListView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'data', 'width'],
      providedOptions: options
    });

    var renderInto = options.renderInto;
    var data = options.data;
    var rowHeight = 30;
    var rowHeightMargin = 5;

    var g = renderInto.append("g")
        .attr("class", "sample-mutation-list");

    g.append("text")
        .text("Chromosome: " + data.chrom);

    data.samples.forEach(function(sample, index) {
      var translateString =
        utils.svgTranslateString(0, (rowHeight + rowHeightMargin) * index);

      var rowContainer = g.append("g")
          .attr("class", "sample-mutation-list__row")
          .attr("transform", translateString);
      
      SampleMutationsView.create({
        renderInto: rowContainer,
        data: data.samples[index],
        width: options.width,
        height: rowHeight,
        chromLength: data.length,
      });
    });
  }

  SampleMutationListView.create = function(options) {
    return new SampleMutationListView(options);
  };

  function SampleMutationsView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'data', 'width', 'height',
        'chromLength'],
      providedOptions: options
    });

    var renderInto = options.renderInto;
    var data = options.data;
    var textWidth = 150;
    var width = options.width - textWidth;
    var height = options.height;
    var chromLength = options.chromLength;

    var xScale = d3.scaleLinear()
      .domain([0, chromLength])
      .range([0, width]);

    var g = renderInto.append("g")
        .attr("class", "sample");

    var label = g.append("text")
        .attr("alignment-baseline", "middle")
        .attr("y", height / 2)
        .text(data.sampleName);

    var background = g.append("rect")
        .attr("transform", utils.svgTranslateString(textWidth, 0))
        .attr("class", "genome-browser__background")
        .attr("width", width)
        .attr("height", height);

    var mutations = g.selectAll(".genome-browser__mutation")
        .data(data.mutationLocations)
      .enter().append("rect")
        .attr("class", "genome-browser__mutation")
        .attr("width", 3)
        .attr("height", height)
        .attr("x", function(d) { return xScale(d) + textWidth; })
        .attr("y", 0);
  }

  SampleMutationsView.create = function(options) {
    return new SampleMutationsView(options);
  };

  return {
    createMutationLocationsView: function(options) {
      return new MutationLocationsView(options);
    }
  };

}(d3, PubSub, utils));
