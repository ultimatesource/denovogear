var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  mutmap.RatioRectangleModel = Backbone.Model.extend({
    defaults: {
      orientation: 'vertical',
      ratioDataList: null
    }
  });

  mutmap.RatioRectangleView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.render();
    },

    render: function() {

      var dim = utils.getD3Dimensions(this.d3el);

      var rectWidth = dim.width;
      var rectHeight = dim.height;

      var data = this.model.get('ratioDataList');

      if (data.length === 0) {
        data = [ { value: 1, color: d3.schemeCategory20[1] } ];
      }

      var vertical = this.model.get('orientation') === 'vertical';

      var partitions = calculatePartitions(data,
        vertical ? dim.height : dim.width);

      var rectUpdate = this.d3el.selectAll('.partition')
          .data(partitions);

      var rectEnter = rectUpdate.enter()
        .append("rect")
          .attr("class", "partition")

      var rectEnterUpdate = rectEnter.merge(rectUpdate);

      rectEnterUpdate
        .attr("x", function(d) {
          return vertical ? 0 : d.offset;
        })
        .attr("y", function(d) {
          return vertical ? d.offset : 0;
        })
        .attr("width", function(d) {
          return vertical ? rectWidth : d.split * rectWidth;
        })
        .attr("height", function(d) {
          return vertical ? d.split * rectHeight : rectHeight;
        })
        .attr("fill", function(d, i) {
          return d.color;
        });
    }
  });

  function calculatePartitions(array, length) {

    var total = _.reduce(array, function(acc, curr) {
      return acc + curr.value;
    }, 0);

    var offset = 0;
    _.each(array, (function(elem) {
      elem.offset = offset;
      elem.split = (elem.value / total);
      offset += (elem.split * length);
    }));

    return array;
  }


}(Backbone, d3, utils));
