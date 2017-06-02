var mutationDistributionView = (function(d3, PubSub, utils) {
  "use strict";

  var optionsManager = utils.createOptionsManager();

  function MutationDistributionView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto'],
      providedOptions: options
    });

    var kinshipPedigreeData = layoutData;
    //var graphData = processPedigree(kinshipPedigreeData);

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

    console.log(width, height);

    svg.append("circle")
        .attr("cx", 0)
        .attr("cy", 0)
        .attr("r", 10)
        .attr("fill", "steelblue");
  }
  
  function createMutationDistributionView(options) {
    return new MutationDistributionView(options);
  }

  return {
    createMutationDistributionView: createMutationDistributionView
  };

}(d3, PubSub, utils));
