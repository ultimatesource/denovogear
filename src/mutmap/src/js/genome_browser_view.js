// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global utils */
/* exported GenomeBrowserView */

var GenomeBrowserView = (function(d3, PubSub) {
  "use strict";

  var GenomeBrowserView = function(selection, vcfData, metadata) {
    this._selection = selection;
    this._vcfData = vcfData;
    this._metadata = metadata;

    this._create();
    this.update();

    PubSub.subscribe("WINDOW_RESIZE", this._onWindowResize.bind(this));
  };

  GenomeBrowserView.prototype._subscribe = function(topic, index) {
    console.log("event received", index);
  };

  GenomeBrowserView.prototype._create = function() {
    this._browser = this._selection.append("svg")
        .attr("class", "genome-browser");

    this._browser.style("width", "100%").style("height", "100%");

    this._margins = {
      left: 20,
      right: 20,
      top: 30,
      bottom: 30
    };

    this._browser.append("rect");

    var translateString = utils.svgTranslateString(this._margins.left,
      this._margins.top);

    this._mutationContainer = this._browser.append("g")
        .attr("transform", translateString)
        .attr("class", "genome-browser__mutation-container")
  }

  GenomeBrowserView.prototype.update = function() {
    var width = parseInt(this._browser.style("width"));
    var height = parseInt(this._browser.style("height"));

    var rectWidth = width - (this._margins.left + this._margins.right);
    this._browser.select("rect")
        .attr("x", this._margins.left)
        .attr("y", this._margins.top)
        .attr("width", rectWidth)
        .attr("height", height - (this._margins.top + this._margins.bottom))
        .style("fill", "#999999");

    var xScale = d3.scaleLinear()
      .domain([this._metadata.minPos, this._metadata.maxPos])
      .range([0, rectWidth]);

    var mutationUpdate = 
      this._mutationContainer.selectAll(".genome-browser__mutation")
        .data(this._vcfData.records);

    var mutationEnter = mutationUpdate.enter().append("rect")
        .attr("class", "genome-browser__mutation")
        .attr("width", 10)
        .attr("height", 100)
        .style("fill", "tomato")
        .on("click", function(d, i) { mutationClicked(d, i); })

    var mutationEnterUpdate = mutationEnter.merge(mutationUpdate);

    mutationEnterUpdate
        .attr("x", function(d) { return xScale(d.POS); })
        .attr("y", 0);
  };

  GenomeBrowserView.prototype._onWindowResize = function() {
    this.update();
  };

  function mutationClicked(d, i) {
    PubSub.publish("MUTATION_CLICKED", { mutationRecordIndex: i });
  };

  return GenomeBrowserView;

}(d3, PubSub));
