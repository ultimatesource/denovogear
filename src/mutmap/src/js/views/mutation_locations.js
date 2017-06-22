
var mutationLocationsView = (function(d3, PubSub, utils) {
  "use strict";

  function MutationLocationsView(options) {

    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'mutationLocationData'],
      providedOptions: options
    });

    var parent = options.renderInto;

    var container = parent.append("div")
        .attr("class", "row")
      .append("div")
        .attr("class",
              "mutation-location-container col-xs-12 panel panel-default")
 
    var chromSelectorContainer = container.append("div")
        .attr("class", "chrom-selector");

    ListSelectorView.create({
      renderInto: chromSelectorContainer,
      list: options.mutationLocationData,
      selector: 'chrom',
      selectedEvent: 'CHROM_SELECTED'
    });

    PubSub.subscribe('CHROM_SELECTED', function(topic, data) {
      console.log(topic, data);
    });

    var svg = container.append("svg")
        .style("width", "100%")
        .style("height", "100%");

    var boundingRect = parent.node().getBoundingClientRect();
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
       SampleMutationsListView.create({
      renderInto: g,
      data: options.mutationLocationData[0],
      width: chromWidth,
    });

  }

  function SampleMutationsListView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'data', 'width'],
      providedOptions: options
    });

    var renderInto = options.renderInto;
    var data = options.data;
    var rowHeight = 30;
    var rowHeightMargin = 5;

    var g = renderInto.append("g")
        .attr("class", "sample-mutation-list");

    this._rows = [];
    data.samples.forEach(function(sample, index) {
      var translateString =
        utils.svgTranslateString(0, (rowHeight + rowHeightMargin) * index);

      var rowContainer = g.append("g")
          .attr("class", "sample-mutation-list__row")
          .attr("transform", translateString);
      
      this._rows.push(SampleMutationsView.create({
        renderInto: rowContainer,
        data: data.samples[index],
        width: options.width,
        height: rowHeight,
        chromLength: data.length,
      }));
    }, this);
  }

  SampleMutationsListView.prototype.update = function(options) {

    optionsManager.checkOptions({
      requiredOptions: ['data'],
      providedOptions: options
    });

    this._rows.forEach(function(sample, index) {
      row.update({
        data: data.samples[index]
      });
    });
  };

  SampleMutationsListView.create = function(options) {
    return new SampleMutationsListView(options);
  };

  function SampleMutationsView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'data', 'width', 'height',
        'chromLength'],
      providedOptions: options
    });

    this._g = options.renderInto.append("g")
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

    this.update(options);
  }

  SampleMutationsView.prototype.update = function(options) {

    optionsManager.checkOptions({
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

    var mutationsEnter = mutationsUpdate.enter();

    mutationsUpdate.exit().remove();
    
    mutationsEnter.append("rect")
        .attr("class", "genome-browser__mutation")
        .attr("width", 3)
        .attr("height", height)
        .attr("y", 0);

    var mutationsEnterUpdate = mutationsEnter.merge(mutationsUpdate);

    mutationsEnterUpdate
        .attr("x", function(d) {
          return xScale(d) + textWidth; 
        });
  };

  SampleMutationsView.create = function(options) {
    return new SampleMutationsView(options);
  };

  function ListSelectorView(options) {
    optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'list', 'selector', 'selectedEvent'],
      providedOptions: options
    });

    var renderInto = options.renderInto;
    var list = options.list;
    var selector = options.selector;

    var currentChromIndex = 0;

    var prevButton = renderInto.append("button")
        .attr("class", "btn btn-default")
        .on("click", function() {
          currentChromIndex--;
          if (currentChromIndex < 0) {
            currentChromIndex = list.length - 1;
          }

          PubSub.publish(options.selectedEvent, currentChromIndex);
        });
    prevButton.append("span")
        .attr("class", "glyphicon glyphicon-arrow-left");
    prevButton.append("span")
        .text(" Previous");

    var dropdown = renderInto.append("div")
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
        .data(list)
      .enter().append("li")
        .attr("class", "mutation")
      .append("a")
        .attr("href", "#")
        .text(function(d) {
          if (selector !== undefined) {
            return d[selector];
          }
          else {
            return d;
          }
        })
        .on("click", function(d, i) {
          currentChromIndex = i;
          PubSub.publish(options.selectedEvent, currentChromIndex);
        });

      var nextButton = renderInto.append("button")
          .attr("class", "btn btn-default")
          .on("click", function() {
            currentChromIndex++;
            if (currentChromIndex === list.length) {
              currentChromIndex = 0;
            }
            PubSub.publish(options.selectedEvent, currentChromIndex);
          });
      nextButton.append("span")
          .text("Next ");
      nextButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-right");


  }

  ListSelectorView.create = function(options) {
    return new ListSelectorView(options);
  };

  return {
    createMutationLocationsView: function(options) {
      return new MutationLocationsView(options);
    }
  };

}(d3, PubSub, utils));
