var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  mutmap.RatioSquareModel = Backbone.Model.extend({
    defaults: {
      ratioDataList: null
    }
  });

  mutmap.RatioSquareView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.render();
    },

    render: function() {

      var dim = utils.getD3Dimensions(this.d3el);

      var squareWidth = dim.width;
      var squareHeight = dim.height;

      var data = this.model.get('ratioDataList');

      if (data.length === 0) {
        data = [ { value: 1, color: d3.schemeCategory20[1] } ];
      }

      var partitions = calculatePartitions(data, dim.width);

      console.log(partitions);

      var squareUpdate = this.d3el.selectAll('.partition')
          .data(partitions);

      var squareEnter = squareUpdate.enter()
        .append("rect")
          .attr("class", "partition")

      var squareEnterUpdate = squareEnter.merge(squareUpdate);

      squareEnterUpdate
        .attr("x", function(d) {
          return d.offset
        })
        .attr("y", 0)
        .attr("width", function(d) {
          return d.split * squareWidth;
        })
        .attr("height", squareHeight)
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
