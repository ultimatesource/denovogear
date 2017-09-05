// eslint exceptions
//
/* exported utils */

var utils = (function() {
  "use strict";

  function halfwayBetween(a, b) {
    if (a > b) {
      return a + ((b - a) / 2);
    }
    else {
      return b + ((a - b) / 2);
    }
  }

  function distanceBetweenPoints(x1, y1, x2, y2) {
    var a = Math.abs(x2 - x1);
    var b = Math.abs(y2 - y1);
    var c = Math.sqrt(Math.pow(a, 2) + Math.pow(b, 2));
    return c;
  }

  function svgTranslateString(x, y) {
    return "translate(" + x + "," + y + ")";
  }

  function sortByKey(options) {

    var array = options.array;
    var key = options.sortKey;
    var descending = options.descending;

    return array.sort(function(a, b) {
      if (a[key] < b[key]) {
        return descending ? 1 : -1;
      }
      if (a[key] > b[key]) {
        return descending ? -1 : 1;
      }
      return 0;
    });
  }

  function getDimensions(element) {

    var node = element.node();
    var boundingRect = node.getBoundingClientRect();

    return {
      width: boundingRect.width,
      height: boundingRect.height
    };
  }

  function getSVGDimensions(element) {
    var bbox = element.node().getBBox();
    return {
      width: bbox.width,
      height: bbox.height
    };
  }

  function getD3Dimensions(element) {
    return {
      width: +element.attr("width"),
      height: +element.attr("height")
    };
  }

  return {
    halfwayBetween: halfwayBetween,
    distanceBetweenPoints: distanceBetweenPoints,
    svgTranslateString: svgTranslateString,
    sortByKey: sortByKey,
    getDimensions: getDimensions,
    getSVGDimensions: getSVGDimensions,
    getD3Dimensions: getD3Dimensions,
    setSizedGroupDimensions: function(width, height) {
      // TODO: Maybe too much of a hack. Forces the g group dimensions by
      // appending an invisible rectangle of the desired size.
      return function(selection) {

        var rectUpdate = selection.selectAll("rect")
          .data([0]);

        var rectEnter = rectUpdate.enter()
          .append("rect")
            .attr("fill", "none");

        var rectEnterUpdate = rectEnter.merge(rectUpdate);

        rectEnterUpdate
            .attr("width", width)
            .attr("height", height);
      };
    },
    getSizedGroupDimensions: function(g) {
      // TODO: same hack as setSizedGroupDimensions
      var rect = g.select("rect");
      return {
        width: +rect.attr("width"),
        height: +rect.attr("height")
      };
    }
  };

}());
