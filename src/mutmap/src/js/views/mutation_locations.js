var mutmap = mutmap || {};

(function(d3, utils) {
  "use strict";

  mutmap.MutationLocationsView = Backbone.View.extend({
    initialize: function(options) {
      this.render(options);
    },

    render: function(options) {

      this.d3el = d3.select(this.el);

      var container = this.d3el.append("div")
          .attr("class", "row")
        .append("div")
          .attr("class",
                "mutation-location-container col-xs-12 panel panel-default")

      var chromSelectorContainer = container.append("div")
          .attr("class", "chrom-selector");

      var SelectedChromosomeModel = Backbone.Model.extend({
        defaults: {
          index: 0
        }
      });

      var selectedChromosome = new SelectedChromosomeModel();

      var chromSelector = new mutmap.ListSelectorView({
        el: chromSelectorContainer,
        model: selectedChromosome,
        list: options.mutationLocationData,
        selector: 'chrom',
        itemName: 'Chromosome'
      });

      selectedChromosome.on('change', function() {
        listView.render({
          data: options.mutationLocationData[selectedChromosome.get('index')]
        });
      });

      var svg = container.append("svg")
          .style("width", "100%")
          .style("height", "100%");

      var boundingRect = this.d3el.node().getBoundingClientRect();
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

      var listView = new SampleMutationsListView({
        el: g,
        data: options.mutationLocationData[0],
        width: chromWidth,
      });

    }
  });

  var SampleMutationsListView = Backbone.View.extend({

    initialize: function(options) {

      this._width = options.width;

      this.render(options);
    },

    render: function(options) {

      utils.optionsManager.checkOptions({
        requiredOptions: ['data'],
        providedOptions: options
      });

      var data = options.data;
      this._rowHeight = 30;
      var rowHeightMargin = 5;

      // TODO: This is a hack. The object-oriented approach I'm taking doesn't
      // really work with d3 all that well. I'm forcing a complete re-render
      // when the data changes. d3 can do more efficient updates, as well as
      // fancy transitions if you stay with d3's model. This works fine for now
      // but might need to be redesigned in the future.
      this.el.selectAll("g").remove();

      var g = this.el.append("g")
          .attr("class", "sample-mutation-list");

      this._rows = [];
      data.samples.forEach(function(sample, index) {
        var translateString =
          utils.svgTranslateString(0,
            (this._rowHeight + rowHeightMargin) * index);

        var rowContainer = g.append("g")
            .attr("class", "sample-mutation-list__row")
            .attr("transform", translateString);
        
        this._rows.push(new SampleMutationsView({
          el: rowContainer,
          data: data.samples[index],
          width: this._width,
          height: this._rowHeight,
          chromLength: data.length,
        }));
      }, this);

      this._rows.forEach(function(sample, index) {
        sample.render({
          data: options.data.samples[index],
          height: this._rowHeight,
          chromLength: options.data.length
        });
      }, this);
    }

  });

  var SampleMutationsView = Backbone.View.extend({

    initialize: function(options) {

      this._g = options.el.append("g")
          .attr("class", "sample");

      var label = this._g.append("text")
          .attr("alignment-baseline", "middle")
          .attr("y", options.height / 2)
          .text(options.data.sampleName);

      this._textWidth = 150;
      this._width = options.width - this._textWidth;

      var background = this._g.append("rect")
          .attr("transform", utils.svgTranslateString(this._textWidth, 0))
          .attr("class", "genome-browser__background")
          .attr("width", this._width)
          .attr("height", options.height);

      this.render(options);
    },

    render: function(options) {

      utils.optionsManager.checkOptions({
        requiredOptions: ['data', 'height',
          'chromLength'],
        providedOptions: options
      });

      var data = options.data;
      var height = options.height;
      var chromLength = options.chromLength;

      var xScale = d3.scaleLinear()
        .domain([0, chromLength])
        .range([0, this._width]);

      // preserve this context
      var textWidth = this._textWidth;

      var mutationsUpdate = this._g.selectAll(".genome-browser__mutation")
          .data(data.mutationLocations);

      var mutationsEnter = mutationsUpdate.enter().append("rect")
          .attr("class", "genome-browser__mutation")
          .attr("width", 3)
          .attr("height", height)
          .attr("y", 0);

      mutationsUpdate.exit().remove();
      

      var mutationsEnterUpdate = mutationsEnter.merge(mutationsUpdate);

      mutationsEnterUpdate
          .attr("x", function(d) {
            return xScale(d) + textWidth; 
          });
    }
  });

}(d3, utils));
