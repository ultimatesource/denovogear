// eslint exceptions
//
/* exported utils */

var utils = (function() {
  "use strict";

  var optionsManager = createOptionsManager();

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

  function OptionsManager() {
  }

  OptionsManager.prototype.checkOptions = function(options) {

    if (options === undefined) throw "No options";
    if (options.requiredOptions  === undefined) throw "No requiredOptions";
    if (options.providedOptions  === undefined) throw "No providedOptions";

    // This is simply used to maintain the list of required options in a more
    // efficient data struct than an actual list, since it will need to be
    // searched for each potential unused option
    //var requiredObject = {};

    // Verify all required options are provided, and also build
    // requiredObject for the unused option check which comes later 
    options.requiredOptions.forEach(function(option) {
      //requiredObject[option] = null;
      if (options.providedOptions[option] === undefined) {
        throw "Required option '" + option + "' not provided";
      }
    });

    // Verify no unused options are provided
    //Object.keys(options.providedOptions).forEach(function(option) {
    //  if (requiredObject[option] === undefined) {
    //    throw "Unused option '" + option + "' provided";
    //  }
    //});
  }

  function createOptionsManager() {
    return new OptionsManager();
  }

  function sortByKey(options) {

    optionsManager.checkOptions({
      requiredOptions: ['array', 'sortKey', 'descending'],
      providedOptions: options
    });

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

  return {
    halfwayBetween: halfwayBetween,
    distanceBetweenPoints: distanceBetweenPoints,
    svgTranslateString: svgTranslateString,
    createOptionsManager: createOptionsManager,
    sortByKey: sortByKey,
    optionsManager: optionsManager
  };

}());
