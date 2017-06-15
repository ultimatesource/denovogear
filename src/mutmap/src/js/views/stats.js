// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* exported statsView */


var statsView = (function(d3, PubSub) {
  "use strict";

  PubSub.subscribe("ACTIVE_NODE_CHANGED", stateChanged);
  PubSub.subscribe("PEDIGREE_VIEW_UPDATE", stateChanged);

  var optionsManager = utils.createOptionsManager();

  var format = d3.format(",.6e");

  var stats = [ "ID", "GT", "GQ", "GP", "DP", "MUP", "MU1P" ];

  var StatsView = function(options) {

    optionsManager.checkOptions({
      requiredOptions: ['renderInto'],
      providedOptions: options
    });

    var stat = options.renderInto.selectAll(".stat")
        .data(stats)
      .enter().append("div")
        .attr("class", "stat");

    stat.append("div")
        .attr("class", "stat__label")
        .text(function(d) { return d + ":"; });

    stat.append("input")
        .attr("class", "stat__text form-control")
        .attr("id", function(d) { return d.toLowerCase() + "_display"; })
        .attr("type", "text")
        .attr("placeholder", function(d) { return d; });
  };

  StatsView.prototype.update = function() {
  };

  function stateChanged(topic, data) {

    if (data.activeNode) {

      var dngData = null;

      var d = data.activeNode;
      if (d.dataNode.data.dngOutputData !== undefined) {
        dngData = d.dataNode.data.dngOutputData;
      }
      else {
        // TODO: should probably be some sort of search for the correct
        // child, rather than assuming it's the first one.
        if (d.dataNode.data.sampleIds.children[0]) {
          dngData = d.dataNode.data.sampleIds.children[0].dngOutputData;
        }
        else {
          dngData = d.dataNode.data.sampleIds.dngOutputData;
        }
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

  function createStatsView(options) {
    return new StatsView(options);
  }

  return {
    createStatsView: createStatsView,
  };

}(d3, PubSub));
