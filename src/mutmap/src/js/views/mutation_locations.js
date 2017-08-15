var mutmap = mutmap || {};

(function(d3, utils) {
  "use strict";


  var SelectedChromosomeModel = Backbone.Model.extend({
    defaults: {
      index: 0
    }
  });

  var SampleListModel = Backbone.Model.extend({
    defaults: {
      listElements: null
    }
  });


  mutmap.MutationLocationsView = Backbone.View.extend({
    initialize: function(options) {

      this.mutationLocationData = options.mutationLocationData;

      this.d3el = d3.select(this.el);

      this.container = this.d3el.append("div");

      this._chromSelectorContainer = this.container.append("div")
          .attr("class", "chrom-selector");

      this._selectedChromosome = new SelectedChromosomeModel();

      this._svg = this.container.append("svg")
          .style("width", "100%");

      this._g = this._svg.append("g");

      var chromSelector = new mutmap.ListSelectorView({
        el: this._chromSelectorContainer,
        model: this._selectedChromosome,
        list: this.mutationLocationData,
        selector: 'chrom',
        itemName: 'Chromosome'
      });

      this._sampleList = new SampleListModel();
      this._sampleList.set('listElements', this.mutationLocationData[0]);

      this._selectedChromosome.on('change', function() {
        this._sampleList.set('listElements',
          this.mutationLocationData[this._selectedChromosome.get('index')]);
      }, this);

      this._g.call(utils.setSizedGroupDimensions(0, 0));

      this._listView = new SampleMutationsListView({
        el: this._g.node(),
        sampleList: this._sampleList,
      });

      this.render();
    },

    render: function() {

      this.dim = utils.getDimensions(this.d3el);

      var selectorDimensions =
        utils.getDimensions(this._chromSelectorContainer);

      var svgHeight = this.dim.height - selectorDimensions.height;
      this._svg.style("height", svgHeight+'px');

      var margins = {
        left: 20,
        top: 20,
        right: 20,
        bottom: 20
      };

      var chromWidth = this.dim.width - (margins.left + margins.right);

      var listWidth = chromWidth;
      var listHeight = svgHeight - (margins.top + margins.bottom);

      this._g.attr("transform",
        utils.svgTranslateString(margins.left, margins.top));

      this._g.call(utils.setSizedGroupDimensions(listWidth, listHeight));

      this._listView.render();

    },

  });


  var SampleMutationsListView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);
      this._sampleList = options.sampleList;

      this._sampleList.on('change', this.render.bind(this));

      this.render();
    },

    render: function() {

      var d3el = d3.select(this.el);
      this.dim = utils.getSizedGroupDimensions(this.d3el);

      var data = this._sampleList.get('listElements');
      var rowHeight = this.dim.height / data.samples.length;
      var rowPadding = 0.15 * rowHeight;

      // TODO: This is a hack. The object-oriented approach I'm taking doesn't
      // really work with d3 all that well. I'm forcing a complete re-render
      // when the data changes. d3 can do more efficient updates, as well as
      // fancy transitions if you stay with d3's model. This works fine for now
      // but might need to be redesigned in the future.
      this.d3el.selectAll("g").remove();

      var g = this.d3el.append("g")
          .attr("class", "sample-mutation-list");

      this._rows = [];
      data.samples.forEach(function(sample, index) {
        var translateString =
          utils.svgTranslateString(0, rowHeight * index);

        var rowContainer = g.append("g")
            .attr("class", "sample-mutation-list__row")
            .attr("transform", translateString);

        var paddedRowHeight = rowHeight - rowPadding;

        rowContainer.call(
          utils.setSizedGroupDimensions(this.dim.width, paddedRowHeight));
        
        this._rows.push(new SampleMutationsView({
          el: rowContainer.node(),
          data: data.samples[index],
          chromLength: data.length,
        }));
      }, this);

      this._rows.forEach(function(sample, index) {
        sample.render();
      }, this);
    }

  });


  var SampleMutationsView = Backbone.View.extend({

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.data = options.data;
      this.chromLength = options.chromLength;

      this._g = this.d3el.append("g")
          .attr("class", "sample");

      this._label = this._g.append("text");

      this._background = this._g.append("rect")
          .attr("class", "genome-browser__background")

      this.render();
    },

    render: function() {

      this.dim = utils.getSVGDimensions(this.d3el);

      var data = this.data;
      var chromLength = this.chromLength;

      var fontSize = Math.min(18, .8 * this.dim.height);
      //this._textWidth = 8 * this.dim.height;
      this._textWidth = 10 * fontSize;
      this._width = this.dim.width - this._textWidth;

      var xScale = d3.scaleLinear()
        .domain([0, chromLength])
        .range([0, this._width]);

      // preserve this context
      var textWidth = this._textWidth;

      this._label
          .attr("alignment-baseline", "middle")
          .attr("font-size", fontSize+"px")
          .attr("y", this.dim.height / 2)
          .text(this.data.sampleName);

      this._background
          .attr("transform", utils.svgTranslateString(this._textWidth, 0))
          .attr("width", this._width)
          .attr("height", this.dim.height);

      var mutationsUpdate = this._g.selectAll(".genome-browser__mutation")
          .data(data.mutationLocations);

      var mutationsEnter = mutationsUpdate.enter().append("rect")
          .attr("class", "genome-browser__mutation")
          .attr("width", 3)
          .attr("height", this.dim.height)
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
