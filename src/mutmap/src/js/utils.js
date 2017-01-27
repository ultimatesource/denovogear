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

  return {
    halfwayBetween: halfwayBetween,
    distanceBetweenPoints: distanceBetweenPoints
  };

}());
