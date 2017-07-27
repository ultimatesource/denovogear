// eslint exceptions
//
/* global d3 */
/* global utils */
/* exported GenomeBrowserView */

var mutmap = mutmap || {};

(function(d3, utils) {
  "use strict";

  mutmap.ContigView = Backbone.View.extend({

    events: {
      'click .prev-mutation-button': 'onPrevMutationButtonClicked',
      'click .next-mutation-button': 'onNextMutationButtonClicked',
      'click .genome-browser__mutation': 'onMutationClicked'
    },

    onPrevMutationButtonClicked: function() {
      this.prevMutationButtonClicked();
    },

    onNextMutationButtonClicked: function() {
      this.nextMutationButtonClicked();
    },

    onMutationClicked: function(e) {

      var selection = d3.select(e.target);
      var d = selection.datum();

      this.mutationClicked(d);
    },

    initialize: function(options) {

      this.d3el = d3.select(this.el);

      this.model.on('change', this.render.bind(this));

      this.prevMutationButtonClicked = options.prevMutationButtonClicked;
      this.nextMutationButtonClicked = options.nextMutationButtonClicked;
      this.mutationClicked = options.mutationClicked;

      this._create();
      this.render();
    },

    _create: function() {

      this._margins = {
        left: 40,
        right: 40,
        top: 5,
        bottom: 60
      };

      var mutationSelector = mutationSelectorMaker();
      this.d3el.call(mutationSelector);

      this._browser = this.d3el.append("svg")
          .attr("class", "genome-browser");

      this._browser.style("width", "100%").style("height", "100%");

      this._rect = this._browser.append("rect")
          .attr("class", "genome-browser__background");

      this._mutationContainer = this._browser.append("g")
          .attr("transform",
            utils.svgTranslateString(this._margins.left, this._margins.top))
          .attr("class", "genome-browser__mutation-container")

      this._mutations = mutationMaker();
    },

    render: function() {

      var boundingRect = this.d3el.node().getBoundingClientRect();
      var width = boundingRect.width;
      var height = boundingRect.height;

      var horizontalMargin = this._margins.left + this._margins.right;
      var verticalMargin = this._margins.top + this._margins.bottom;

      var genomeWidth = width - horizontalMargin;
      var genomeHeight = height - verticalMargin;

      this._rect
          .attr("transform",
            utils.svgTranslateString(this._margins.left, this._margins.top))
          .attr("x", 0)
          .attr("y", 0)
          .attr("width", genomeWidth)
          .attr("height", genomeHeight);

      var contigLength = this.model.get('contigData').length;
      this._xFocusScale = d3.scaleLinear()
        .domain([0, contigLength])
        .range([0, genomeWidth]);
      this._xContextScale = d3.scaleLinear()
        .domain(this._xFocusScale.domain())
        .range(this._xFocusScale.range());

      this._xAxis = d3.axisBottom(this._xFocusScale);

      // TODO: Hack. Find a proper way to update the original axes without
      // forcibly deleting the old ones
      this._browser.selectAll(".axis--x").remove();
      this._browser.append("g")
          .attr("class", "axis axis--x")
          .attr("transform", utils.svgTranslateString(
            this._margins.left, genomeHeight + this._margins.top))
          .call(this._xAxis);

      //this._mutationContainer.selectAll(".genome-browser__mutation").remove();
      this._mutations
        .vcfData(this.model.get('contigData'))
        .scale(this._xFocusScale)
        .height(genomeHeight)
        .selectedMutationIndex(this.model.get('selectedMutationIndex'));
      this._mutationContainer.call(this._mutations);

    },

    _onWindowResize: function() {
      this.render();
    },

  });

  function mutationMaker() {

    var scale;
    var vcfData;
    var height;
    var selectedMutationIndex;

    function my(selection) {

      var mutationWidth = 6;

      var mutationUpdate = 
        selection.selectAll(".genome-browser__mutation")
          .data(vcfData.records);

      var mutationEnter = mutationUpdate.enter().append("rect")
          .attr("width", mutationWidth)
          .attr("height", height);

      var mutationExit = mutationUpdate.exit();
      mutationExit.remove();

      var mutationEnterUpdate = mutationEnter.merge(mutationUpdate);

      // preserve 'this' context
      var xFocusScale = scale;
      mutationEnterUpdate
          .attr("class", function(d, i) {
            if (i === selectedMutationIndex) {
              return "genome-browser__mutation genome-browser__mutation--selected";
            }
            else {
              return "genome-browser__mutation";
            }
          })
          .attr("x", function(d) { return xFocusScale(d.POS); })
          .attr("y", 0);

      mutationEnterUpdate.each(function(d, i) {
        if (i === selectedMutationIndex) {
          // Make this the last element in the DOM parent so it renders on top
          // of any overlapping mutations
          d3.select(this).raise();
        }
      });
    }

    my.scale = function(value) {
      if (!arguments.length) return scale;
      scale = value;
      return my;
    };

    my.vcfData = function(value) {
      if (!arguments.length) return vcfData;
      vcfData = value;
      return my;
    };

    my.height = function(value) {
      if (!arguments.length) return height;
      height = value;
      return my;
    };

    my.selectedMutationIndex = function(value) {
      if (!arguments.length) return selectedMutationIndex ;
      selectedMutationIndex = value;
      return my;
    };

    return my;
  }

  function mutationSelectorMaker() {

    function my(selection) {

      var prevButton = selection.append("button")
          .attr("class", "prev-mutation-button btn btn-default")
      prevButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-left");
      prevButton.append("span")
          .text(" Previous Mutation");

      var nextButton = selection.append("button")
          .attr("class", "next-mutation-button btn btn-default")
      nextButton.append("span")
          .text("Next Mutation ");
      nextButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-right");
    }

    return my;
  }

}(d3, utils));
