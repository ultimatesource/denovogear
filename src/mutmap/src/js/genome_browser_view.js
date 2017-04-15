// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global utils */
/* exported GenomeBrowserView */

var GenomeBrowserView = (function(d3, PubSub) {
  "use strict";

  var GenomeBrowserView = function(selection, vcfData) {
    this._selection = selection;
    this._vcfData = vcfData;

    this._create();
    this.update();

    PubSub.subscribe("WINDOW_RESIZE", this._onWindowResize.bind(this));
  };

  GenomeBrowserView.prototype._create = function() {
    this._browser = this._selection.append("svg")
        .attr("class", "genome-browser");

    this._browser.style("width", "100%").style("height", "100%");

    this._margins = {
      left: 40,
      right: 40,
      top: 5,
      bottom: 20
    };

    this._focus = this._browser.append("g")
          .attr("class", "genome-browser__focus")
      
    this._focusRect = this._focus.append("rect")
        .attr("class", "genome-browser__background");

    this._context = this._browser.append("g")
          .attr("class", "genome-browser__context")
      
    this._contextRect = this._context.append("rect")
        .attr("class", "genome-browser__brush");

    var translateString = utils.svgTranslateString(this._margins.left,
      this._margins.top);

    this._mutationContainer = this._browser.append("g")
        .attr("transform", translateString)
        .attr("class", "genome-browser__mutation-container")

    this._contextMutationContainer = this._browser.append("g");
  }

  GenomeBrowserView.prototype.update = function() {
    var width = parseInt(this._browser.style("width"));
    var height = parseInt(this._browser.style("height"));

    var verticalMargin = this._margins.top + this._margins.bottom;
    var marginedHeight = height - verticalMargin; 
    var focusHeight = 0.5 * marginedHeight;
    var xAxisPadding = 20;
    var contextHeight = marginedHeight - focusHeight - xAxisPadding;

    this._focusRect
        .attr("transform",
          utils.svgTranslateString(this._margins.left, this._margins.top));

    var rectWidth = width - (this._margins.left + this._margins.right);
    this._focusRect
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", rectWidth)
        .attr("height", focusHeight);

    var contigLength = this._vcfData.header.contig[0]['length'];
    this._xFocusScale = d3.scaleLinear()
      .domain([0, contigLength])
      .range([0, rectWidth]);
    this._xContextScale = d3.scaleLinear()
      .domain(this._xFocusScale.domain())
      .range(this._xFocusScale.range());

    this._xAxis = d3.axisBottom(this._xFocusScale);
    var xAxisContext = d3.axisBottom(this._xContextScale);

    var mutations = mutationMaker()
      .vcfData(this._vcfData)
      .scale(this._xFocusScale)
      .height(focusHeight);
    this._mutationContainer.call(mutations);

    var contextYOffset = this._margins.top + focusHeight + xAxisPadding;

    this._contextMutationContainer
        .attr("transform", utils.svgTranslateString(
          this._margins.left, contextYOffset))
        .attr("class", "genome-browser__mutation-container")

    this._contextMutationContainer.call(mutations.height(contextHeight));

    this._context
        .attr("transform", utils.svgTranslateString(this._margins.left,
          contextYOffset));

    this._contextRect
        .attr("class", "genome-browser__background")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", rectWidth)
        .attr("height", contextHeight);

    var brush = d3.brushX()
      .extent([[0, 0], [rectWidth, contextHeight]])
      .on("brush end", this._brushed.bind(this));

    // TODO: Hack. Find a proper way to update the original brush without
    // forcibly deleting the old one
    this._context.select(".brush").remove();

    this._context.append("g")
        .attr("class", "brush")
        .call(brush)
        // Start brush covering entire genome
        .call(brush.move, this._xFocusScale.range())

    // TODO: Hack. Find a proper way to update the original axes without
    // forcibly deleting the old ones
    this._focus.selectAll(".axis--x").remove();

    this._focus.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", utils.svgTranslateString(
          this._margins.left, focusHeight + this._margins.top))
        .call(this._xAxis);

    this._focus.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", utils.svgTranslateString(
          this._margins.left, focusHeight + this._margins.top +
            xAxisPadding + contextHeight))
        .call(this._xAxis);
  };

  GenomeBrowserView.prototype._brushed = function() {
    if (d3.event.sourceEvent && d3.event.sourceEvent.type === "zoom") return;

      var s = d3.event.selection || this._xContextScale.range();
      this._xFocusScale.domain(s.map(this._xContextScale.invert,
        this._xContextScale));

      this._focus.select(".axis--x").call(this._xAxis);

      var xFocusScale = this._xFocusScale;
      this._mutationContainer.selectAll(".genome-browser__mutation")
          .attr("x", function(d) {
            var pos = d.POS;
            var domain = xFocusScale.domain();
            if (pos >= domain[0] && pos <= domain[1]) {
              return xFocusScale(d.POS);
            }
            else {
              // TODO: Hack. Should probably make this invisibly instead of
              // just moving it way off the screen
              return -1000;
            }
          });
  }

  GenomeBrowserView.prototype._onWindowResize = function() {
    this.update();
  };

  function mutationMaker(scale, vcfData) {

    var scale;
    var vcfData;
    var height;

    function my(selection) {
      var mutationWidth = 6;

      var mutationUpdate = 
        selection.selectAll(".genome-browser__mutation")
          .data(vcfData.records);

      var mutationEnter = mutationUpdate.enter().append("rect")
          .attr("class", "genome-browser__mutation")
          .attr("width", mutationWidth)
          .attr("height", height)
          .on("click", function(d, i) { mutationClicked(d, i); })

      var mutationEnterUpdate = mutationEnter.merge(mutationUpdate);

      // preserve 'this' context
      var xFocusScale = scale;
      mutationEnterUpdate
          .attr("x", function(d) { return xFocusScale(d.POS); })
          .attr("y", 0);
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

    return my;
  }

  function mutationClicked(d, i) {
    PubSub.publish("MUTATION_CLICKED", { mutationRecordIndex: i });
  };

  return GenomeBrowserView;

}(d3, PubSub));
