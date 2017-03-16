var genomeBrowserView = (function($, d3) {
  "use strict";

  var vcfData;
  var metadata;

  function createGenomeBrowser() {

    function my(selection) {

      var contigLength = vcfData.header.contig[0].length;

      var svg = selection.append("svg")
          .attr("class", "browser-svg");
      var width = $(".browser-svg").parent().width();
      var height = $(".browser-svg").parent().height();

      svg.attr("width", width).attr("height", height);
      
      var browser = svg.append("g")
          .attr("class", "genomeBrowser");

      var margins = {
        left: 50,
        right: 50,
        top: 40,
        bottom: 40
      };

      var rectWidth = width - (margins.left + margins.right);
      browser.append("rect")
          .attr("x", margins.left)
          .attr("y", margins.top)
          .attr("width", rectWidth)
          .attr("height", height - (margins.top + margins.bottom))
          .style("fill", "#999999");

      var translateString = utils.svgTranslateString(margins.left,
        margins.top);

      //var xScale = d3.scaleLinear()
      //  .domain([0, contigLength])
      //  .range([0, rectWidth]);

      var xScale = d3.scaleLinear()
        .domain([metadata.minPos, metadata.maxPos])
        .range([0, rectWidth]);

      var mutations = browser.append("g")
          .attr("transform", translateString)
          .selectAll(".mutation")
          .data(vcfData.records)
          .enter()
          .append("rect")
          .attr("class", "mutation")
          .attr("x", function(d) { return xScale(d.POS); })
          .attr("y", 0)
          .attr("width", 10)
          .attr("height", 100)
          .style("fill", "tomato")
          .on("click", function(d, i) { mutationClicked(d, i); });
    }

    my.vcfData = function(value) {
      if (!arguments.length) return vcfData;
      vcfData = value;
      return my;
    };

    my.metadata = function(value) {
      if (!arguments.length) return vcfData;
      metadata = value;
      return my;
    };

    function mutationClicked(d, i) {
      var action = {
        type: "MUTATION_CLICKED",
        mutationRecord: d,
        mutationRecordIndex: i
      };
      store.dispatch(action);
    }

    return my;
  }

  return {
    createGenomeBrowser: createGenomeBrowser
  };

}($, d3, store));
