
var mutationLocationsView = (function(d3, PubSub, utils) {
  "use strict";

  function MutationLocationsView(options) {

    var optionsManager = utils.createOptionsManager();

    optionsManager.checkOptions({
      requiredOptions: ['renderInto'],
      providedOptions: options
    });

    var parent = options.renderInto;

    var g = parent.append("div")
        .attr("class", "panel panel-default")
      .append("svg").append("g");
    
    g.append("circle")
        .attr("r", 20)
        .attr("cx", 0)
        .attr("cy", 0);

  }

  return {
    createMutationLocationsView: function(options) {
      return new MutationLocationsView(options);
    }
  };

}(d3, PubSub, utils));
