var mutmap = mutmap || {};

(function(Backbone, d3, utils) {
  "use strict";

  mutmap.RatioCircleModel = Backbone.Model.extend({
    defaults: {
      ratioDataList: null
    }
  });

  mutmap.RatioCircleView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.container = this.d3el.append("g")
          .attr("class", "ratio-circle-container");

      this.render();
    },

    render: function() {

      var dim = utils.getD3Dimensions(this.d3el);

      var data = this.model.get('ratioDataList');

      if (data.length === 0) {
        data = [ { value: 1, color: d3.schemeCategory20[7] } ];
      }

      this.container
          .attr("transform",
            utils.svgTranslateString(dim.width/2, dim.height/2));

      var pie = d3.pie()
        .sort(null)
        .value(function(d) {
          return d.value;
        });
      var piePath = d3.arc()
        .outerRadius(dim.width / 2)
        .innerRadius(0);

      var pieUpdate = this.container.selectAll(".arc")
          .data(pie(data));

      var pieEnter = pieUpdate.enter().append("path")
          .attr("class", "arc");

      var pieEnterUpdate = pieEnter.merge(pieUpdate);

      pieEnterUpdate
          .attr("d", piePath)
          .style("fill", function(d) { 
            return d.data.color;
          });
    }
  });

}(Backbone, d3, utils));
