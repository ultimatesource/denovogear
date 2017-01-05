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

  return { halfwayBetween: halfwayBetween };

}());
