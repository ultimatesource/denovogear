// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global utils */
/* exported GenomeBrowserView */

var mutmap = mutmap || {};

(function(d3, PubSub, utils) {
  "use strict";

  mutmap.ContigView = Backbone.View.extend({

    initialize: function(options) {
      this._vcfData = options.vcfData;
      this._contigData = options.contigData;
      this._selectedContigIndex = options.selectedContigIndex;
      this._selectedMutationIndex = options.selectedMutationIndex;

      this.model.on('change', function() {
      });

      PubSub.subscribe("WINDOW_RESIZE", this._onWindowResize.bind(this));
      PubSub.subscribe("MUTATION_INDEX_UPDATED",
        this._selectedMutationIndexChanged.bind(this));

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

      var mutationSelector = mutationSelectorMaker().vcfData(this._vcfData);
      this.el.call(mutationSelector);

      this._browser = this.el.append("svg")
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

      var boundingRect = this.el.node().getBoundingClientRect();
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

      var contigLength = this._contigData[this._selectedContigIndex].length;
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
        .vcfData(this._contigData[this._selectedContigIndex])
        .scale(this._xFocusScale)
        .height(genomeHeight)
        .selectedMutationIndex(this._selectedMutationIndex);
      this._mutationContainer.call(this._mutations);

    },

    _onWindowResize: function() {
      this.render();
    },

    _selectedMutationIndexChanged: function(topic, data) {

      this._selectedContigIndex = data.selectedContigIndex;
      this._selectedMutationIndex = data.selectedMutationIndex;
      this.render();
    }

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
          .attr("height", height)
          .on("click", mutationClicked)

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

    var vcfData;

    function my(selection) {

      var prevButton = selection.append("button")
          .attr("class", "btn btn-default")
          .on("click", function() {
            PubSub.publish("PREV_MUTATION_BUTTON_CLICKED");
          });
      prevButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-left");
      prevButton.append("span")
          .text(" Previous Mutation");

      var dropdown = selection.append("div")
          .attr("class", "dropdown")
      dropdown.append("button")
          .attr("class", "btn btn-default dropdown-toggle")
          .attr("id", "dropdownMenu1")
          .attr("type", "button")
          .attr("data-toggle", "dropdown")
          .text("Select");

      dropdown.append("ul")
          .attr("class", "dropdown-menu")
        .selectAll(".mutation")
          .data(vcfData.records)
        .enter().append("li")
          .attr("class", "mutation")
        .append("a")
          .attr("href", "#")
          .text(function(d) {
            return "Contig: " + d.CHROM + ", Position: " + d.POS;
          })
          .on("click", function(d) {
            PubSub.publish("MUTATION_SELECTED", d);
          });

      var nextButton = selection.append("button")
          .attr("class", "btn btn-default")
          .on("click", function() {
            PubSub.publish("NEXT_MUTATION_BUTTON_CLICKED");
          });
      nextButton.append("span")
          .text("Next Mutation ");
      nextButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-right");
    }

    my.vcfData = function(value) {
      if (!arguments.length) return vcfData;
      vcfData = value;
      return my;
    };

    return my;
  }

  function mutationClicked(d, i) {
    PubSub.publish("MUTATION_SELECTED", d);
  };

}(d3, PubSub, utils));
