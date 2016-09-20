(function() {
  "use strict";

  var vcf = require('vcf.js');

  function parseVCFText(text) {
    var parsed = vcf.parser()(text);
    var inheritancePattern = parsed.records[0].INFO.DNT;
    var parent0 = inheritancePattern.slice(0, 2);
    var parent1 = inheritancePattern.slice(3, 5);
    var child = inheritancePattern.slice(6, 8);
  }

  module.exports.parseVCFText = parseVCFText;

}());
