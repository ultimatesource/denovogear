// eslint exceptions
//
/* global Backbone */
/* global d3 */

var mutmap = mutmap || {};

(function(Backbone, d3) {
  "use strict";

  var format = d3.format(",.6e");

  var stats = [ "ID", "GT", "GQ", "GP", "DP", "MUP", "MU1P" ];

  mutmap.StatsView = Backbone.View.extend({

    initialize: function() {

      this.model.on('change', this.render.bind(this));

      this.render();
    },

    render: function() {

      var statsUpdate = this.el.selectAll(".stat")
          .data(stats);

      var statsEnter = statsUpdate.enter().append("div")
          .attr("class", "stat");

      statsEnter.append("div")
          .attr("class", "stat__label")
          .text(function(d) { return d + ":"; });

      statsEnter.append("input")
          .attr("class", function(d) {
            return "stat__text form-control " + d.toLowerCase() + "_display"
          })
          .attr("type", "text")
          .attr("placeholder", function(d) { return d; });

      var statsEnterUpdate = statsEnter.merge(statsUpdate);

      this._textBoxes = statsEnterUpdate.selectAll('.stat__text');
     
      var dataNode = this.model.get('dataNode');

      if (dataNode) {

        var dngData = dngDataFromDataNode(dataNode);

        this.el.select(".id_display").attr("value", dataNode.id);
        this.el.select(".gt_display").attr("value", dngData.GT);
        this.el.select(".gq_display").attr("value", dngData.GQ);
        this.el.select(".gp_display").attr("value", dngData.GP);
        this.el.select(".dp_display").attr("value", dngData.DP);
        this.el.select(".mup_display").attr("value", format(dngData.MUP));
        this.el.select(".mu1p_display").attr("value", format(dngData.MU1P));
      }
    },
  });

  function dngDataFromDataNode(dataNode) {

    var dngData = null;

    if (dataNode.data.dngOutputData !== undefined) {
      dngData = dataNode.data.dngOutputData;
    }
    else {
      // TODO: should probably be some sort of search for the correct
      // child, rather than assuming it's the first one.
      if (dataNode.data.sampleIds.children[0]) {
        dngData = dataNode.data.sampleIds.children[0].dngOutputData;
      }
      else {
        dngData = dataNode.data.sampleIds.dngOutputData;
      }
    }

    return dngData;
  }

}(Backbone, d3));
