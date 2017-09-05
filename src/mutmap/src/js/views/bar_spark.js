var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  mutmap.BarSparkModel = Backbone.Model.extend({
    defaults: {
      barDataList: null
    }
  });

  mutmap.BarSparkView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.render();
    },

    render: function() {

      var dim = utils.getD3Dimensions(this.d3el);

      var barDataList = this.model.get('barDataList');

      var max = _.max(barDataList, function(x) { return x.value; }).value;

      var xScale = d3.scaleLinear()
        .domain([0, barDataList.length])
        .range([0, dim.width]);

      var yScale = d3.scaleLinear()
        .domain([0, max])
        .range([dim.height, 0]);

      var barsUpdate = this.d3el.selectAll(".bar")
          .data(this.model.get('barDataList'))

      var barsEnter = barsUpdate.enter();
      var padding = 0.05 * dim.width;

      barsEnter.append("rect")
          .attr("class", "bar")
          .attr("x", function(d, i) {
            return xScale(i);
          })
          .attr("y", function(d, i) {
            return yScale(d.value);
          })
          .attr("width", (dim.width / barDataList.length) - padding)
          .attr("height", function(d) {
            return dim.height - yScale(d.value);
          })
          .style("fill", function(d) {
            return d.color;
          })

    }
  })

}(Backbone, d3, utils));
