// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global utils */
/* exported GenomeBrowserView */

var contigView = (function(d3, PubSub, utils) {
  "use strict";

  var optionsManager = utils.createOptionsManager();

  function ContigView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'vcfData', 'contigData', 
        'selectedContigIndex', 'selectedMutationIndex'],
      providedOptions: options
    });

    this._selection = options.renderInto;
    this._vcfData = options.vcfData;
    this._contigData = options.contigData;
    this._selectedContigIndex = options.selectedContigIndex;
    this._selectedMutationIndex = options.selectedMutationIndex;

    PubSub.subscribe("WINDOW_RESIZE", this._onWindowResize.bind(this));
    PubSub.subscribe("MUTATION_INDEX_UPDATED",
      this._selectedMutationIndexChanged.bind(this));
  }

  var SimpleContigView = function(options) {
    ContigView.call(this, options);
    this._create();
    this.update();
  };

  SimpleContigView.prototype = Object.create(ContigView.prototype);
  SimpleContigView.prototype.constructor = SimpleContigView;

  SimpleContigView.prototype._create = function() {

    this._margins = {
      left: 40,
      right: 40,
      top: 5,
      bottom: 60
    };

    var mutationSelector = mutationSelectorMaker().vcfData(this._vcfData);
    this._selection.call(mutationSelector);

    this._browser = this._selection.append("svg")
        .attr("class", "genome-browser");

    this._browser.style("width", "100%").style("height", "100%");

    this._rect = this._browser.append("rect")
        .attr("class", "genome-browser__background");

    this._mutationContainer = this._browser.append("g")
        .attr("transform",
          utils.svgTranslateString(this._margins.left, this._margins.top))
        .attr("class", "genome-browser__mutation-container")

    this._mutations = mutationMaker();
  };

  SimpleContigView.prototype.update = function() {

    var width = parseInt(this._browser.style("width"));
    var height = parseInt(this._browser.style("height"));

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

  };


  //var GenomeBrowserView = function(options) {
  //  ContigView.call(this, options);

  //  this._create();
  //  this.update();

  //};
  //GenomeBrowserView.prototype = Object.create(ContigView.prototype);
  //GenomeBrowserView.prototype.constructor = GenomeBrowserView;

  //GenomeBrowserView.prototype._create = function() {
  //  this._browser = this._selection.append("svg")
  //      .attr("class", "genome-browser");

  //  this._browser.style("width", "100%").style("height", "100%");

  //  this._margins = {
  //    left: 40,
  //    right: 40,
  //    top: 5,
  //    bottom: 20
  //  };

  //  this._focus = this._browser.append("g")
  //        .attr("class", "genome-browser__focus")
  //    
  //  this._focusRect = this._focus.append("rect")
  //      .attr("class", "genome-browser__background");

  //  this._context = this._browser.append("g")
  //        .attr("class", "genome-browser__context")
  //    
  //  this._contextRect = this._context.append("rect")
  //      .attr("class", "genome-browser__brush");

  //  var translateString = utils.svgTranslateString(this._margins.left,
  //    this._margins.top);

  //  // TODO: Append the mutations to this._focus instead of this._browser,
  //  // in order to match the hierarchical model of the system. I'm leaving it
  //  // this way for now since doing it the other way breaks click events for
  //  // the mutations because of the zooming and brushing and I don't want to
  //  // spend the time to figure that out at the moment.
  //  this._mutationContainer = this._browser.append("g")
  //      .attr("transform", translateString)
  //      .attr("class", "genome-browser__mutation-container")

  //  this._contextMutationContainer = this._browser.append("g");
  //}

  //GenomeBrowserView.prototype.update = function() {
  //  var width = parseInt(this._browser.style("width"));
  //  var height = parseInt(this._browser.style("height"));

  //  var verticalMargin = this._margins.top + this._margins.bottom;
  //  var marginedHeight = height - verticalMargin; 
  //  var focusHeight = 0.5 * marginedHeight;
  //  var xAxisPadding = 20;
  //  var contextHeight = marginedHeight - focusHeight - xAxisPadding;

  //  this._focusRect
  //      .attr("transform",
  //        utils.svgTranslateString(this._margins.left, this._margins.top));

  //  this._genomeWidth = width - (this._margins.left + this._margins.right);
  //  this._focusRect
  //      .attr("x", 0)
  //      .attr("y", 0)
  //      .attr("width", this._genomeWidth)
  //      .attr("height", focusHeight);

  //  var contigLength = this._vcfData.header.contig[0]['length'];
  //  this._xFocusScale = d3.scaleLinear()
  //    .domain([0, contigLength])
  //    .range([0, this._genomeWidth]);
  //  this._xContextScale = d3.scaleLinear()
  //    .domain(this._xFocusScale.domain())
  //    .range(this._xFocusScale.range());

  //  this._xAxis = d3.axisBottom(this._xFocusScale);
  //  var xAxisContext = d3.axisBottom(this._xContextScale);

  //  var mutations = mutationMaker()
  //    .vcfData(this._vcfData)
  //    .scale(this._xFocusScale)
  //    .height(focusHeight);
  //  this._mutationContainer.call(mutations);

  //  var contextYOffset = this._margins.top + focusHeight + xAxisPadding;

  //  this._contextMutationContainer
  //      .attr("transform", utils.svgTranslateString(
  //        this._margins.left, contextYOffset))
  //      .attr("class", "genome-browser__mutation-container")

  //  this._contextMutationContainer.call(mutations.height(contextHeight));

  //  this._context
  //      .attr("transform", utils.svgTranslateString(this._margins.left,
  //        contextYOffset));

  //  this._contextRect
  //      .attr("class", "genome-browser__background")
  //      .attr("x", 0)
  //      .attr("y", 0)
  //      .attr("width", this._genomeWidth)
  //      .attr("height", contextHeight);

  //  this._brush = d3.brushX()
  //    .extent([[0, 0], [this._genomeWidth, contextHeight]])
  //    .on("brush end", this._brushed.bind(this));
  //  
  //  this._zoom = d3.zoom()
  //    .scaleExtent([1, Infinity])
  //    .translateExtent([[0, 0], [this._genomeWidth, focusHeight]])
  //    .extent([[0, 0], [this._genomeWidth, focusHeight]])
  //    .on("zoom", this._zoomed.bind(this));

  //  // TODO: Hack. Find a proper way to update the original brush without
  //  // forcibly deleting the old one
  //  this._context.select(".brush").remove();
  //  this._context.append("g")
  //      .attr("class", "brush")
  //      .call(this._brush)
  //      // Start brush covering entire genome
  //      .call(this._brush.move, this._xFocusScale.range())

  //  // TODO: Hack. Find a proper way to update the original axes without
  //  // forcibly deleting the old ones
  //  this._focus.selectAll(".axis--x").remove();
  //  this._focus.append("g")
  //      .attr("class", "axis axis--x")
  //      .attr("transform", utils.svgTranslateString(
  //        this._margins.left, focusHeight + this._margins.top))
  //      .call(this._xAxis);

  //  this._context.selectAll(".axis--x").remove();
  //  this._context.append("g")
  //      .attr("class", "axis axis--x")
  //      .attr("transform", utils.svgTranslateString(0, contextHeight))
  //      .call(this._xAxis);

  //  this._focus.selectAll(".zoom").remove();
  //  this._focus.append("rect")
  //      .attr("class", "zoom")
  //      .attr("width", this._genomeWidth)
  //      .attr("height", focusHeight)
  //      .attr("transform", "translate(" + this._margins.left + "," +
  //        this._margins.top + ")")
  //      .call(this._zoom);
  //};

  //GenomeBrowserView.prototype._brushed = function() {
  //  // ignore brush-by-zoom
  //  if (d3.event.sourceEvent && d3.event.sourceEvent.type === "zoom") return;

  //  var s = d3.event.selection || this._xContextScale.range();
  //  this._xFocusScale.domain(s.map(this._xContextScale.invert,
  //    this._xContextScale));

  //  this._focus.select(".axis--x").call(this._xAxis);

  //  this._browser.select(".zoom").call(this._zoom.transform, d3.zoomIdentity
  //    .scale(this._genomeWidth / (s[1] - s[0]))
  //    .translate(-s[0], 0));

  //  this._repositionMutations();

  //}

  //GenomeBrowserView.prototype._zoomed = function() {
  //  // ignore zoom-by-brush
  //  if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return;

  //  var t = d3.event.transform;
  //  this._xFocusScale.domain(t.rescaleX(this._xContextScale).domain());
  //  this._focus.select(".axis--x").call(this._xAxis);
  //  this._context.select(".brush").call(this._brush.move,
  //    this._xFocusScale.range().map(t.invertX, t));

  //  this._repositionMutations();

  //};

  //GenomeBrowserView.prototype._repositionMutations = function() {
  //  var xFocusScale = this._xFocusScale;
  //  this._mutationContainer.selectAll(".genome-browser__mutation")
  //      .attr("x", function(d) {
  //        var pos = d.POS;
  //        var domain = xFocusScale.domain();
  //        if (pos >= domain[0] && pos <= domain[1]) {
  //          return xFocusScale(d.POS);
  //        }
  //        else {
  //          // TODO: Hack. Should probably make this invisibly instead of
  //          // just moving it way off the screen
  //          return -1000;
  //        }
  //      });
  //};

  ContigView.prototype._onWindowResize = function() {
    this.update();
  };

  ContigView.prototype._selectedMutationIndexChanged =
      function(topic, data) {
    this._selectedContigIndex = data.selectedContigIndex;
    this._selectedMutationIndex = data.selectedMutationIndex;
    this.update();
  };

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

  function createContigView(options) {
    return new SimpleContigView(options);
    //return new GenomeBrowserView(options);
  }

  return {
    createContigView: createContigView
  };

}(d3, PubSub, utils));
