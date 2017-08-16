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

  var SampleMutationsModel = Backbone.Model.extend({
    defaults: {
      chromLength: 0,
      sample: null
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

      this._g = this.d3el.append("g")
          .attr("class", "sample-mutation-list");

      this._rowViews = [];

      this.render();
    },

    render: function() {

      var d3el = d3.select(this.el);
      this.dim = utils.getSizedGroupDimensions(this.d3el);

      var data = this._sampleList.get('listElements');
      var rowHeight = this.dim.height / data.samples.length;
      var rowPadding = 0.15 * rowHeight;
      var paddedRowHeight = rowHeight - rowPadding;

      var rowUpdate = this._g.selectAll(".sample-mutation-list__row")
          .data(data.samples);

      var rowEnter = rowUpdate.enter()
        .append("g")
          .attr("class", "sample-mutation-list__row");

      // Basically the approach I'm taking combines the flow of d3 with the
      // object oriented view architecture I'm using with Backbone. I'm trying
      // to stick with the enter, update, exit paradigm from d3. But I also
      // need to manage the lifetimes of the child Backbone views. So I'm
      // using this._rowViews to hold those references but using the d3
      // selections to manage them.

      var self = this;

      rowEnter.each(function(row, index) {

        var rowContainer = d3.select(this);

        var model = new SampleMutationsModel({
          sample: data.samples[index],
          chromLength: data.length
        });

        self._rowViews[index] = new SampleMutationsView({
          el: rowContainer.node(),
          model: model,
        });
      })

      rowUpdate.each(function(row, index) {

        self._rowViews[index].model.set({
          sample: data.samples[index],
          chromLength: data.length
        });
      })

      var rowEnterUpdate = rowEnter.merge(rowUpdate);

      rowEnterUpdate
          .attr("transform", function(d, i) {
            return utils.svgTranslateString(0, rowHeight * i);
          })
          .call(
            utils.setSizedGroupDimensions(this.dim.width, paddedRowHeight))
          .each(function(row, index) {
            self._rowViews[index].render();
          })

      var rowExit = rowUpdate.exit();
      rowExit
      //rowExit.each(function(row, index) {
      //  self._rowViews.splice(index, 1);
      //})
        .remove()

    }
  });


  var SampleMutationsView = Backbone.View.extend({

    initialize: function() {

      this.d3el = d3.select(this.el);

      this._g = this.d3el.append("g")
          .attr("class", "sample");

      this._label = this._g.append("text");

      this._background = this._g.append("rect")
          .attr("class", "genome-browser__background")

      this.model.on('change', this.render.bind(this));

      this.render();
    },

    render: function() {

      console.log("render");
      this.dim = utils.getSizedGroupDimensions(this.d3el);

      var fontSize = Math.min(18, .8 * this.dim.height);
      this._textWidth = 10 * fontSize;
      this._width = this.dim.width - this._textWidth;

      var xScale = d3.scaleLinear()
        .domain([0, this.model.get('chromLength')])
        .range([0, this._width]);

      // preserve this context
      var textWidth = this._textWidth;

      this._label
          .attr("alignment-baseline", "middle")
          .attr("font-size", fontSize+"px")
          .attr("y", this.dim.height / 2)
          .text(this.model.get('sample').sampleName);

      this._background
          .attr("transform", utils.svgTranslateString(this._textWidth, 0))
          .attr("width", this._width)
          .attr("height", this.dim.height);

      var mutationsUpdate = this._g.selectAll(".genome-browser__mutation")
          .data(this.model.get('sample').mutationLocations);

      var mutationsEnter = mutationsUpdate.enter().append("rect")
          .attr("class", "genome-browser__mutation")
          .attr("y", 0);

      mutationsUpdate.exit().remove();
      
      var mutationsEnterUpdate = mutationsEnter.merge(mutationsUpdate);

      mutationsEnterUpdate
          .attr("x", function(d) {
            return xScale(d) + textWidth; 
          })
          .attr("width", 3)
          .attr("height", this.dim.height)
    }
  });

}(d3, utils));
