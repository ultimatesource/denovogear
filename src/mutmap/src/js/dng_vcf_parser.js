(function() {
  "use strict";

  var vcf = require('vcf.js');

  function parseVCFText(text) {
    var parsed = vcf.parser()(text);
    return parsed;
    //console.log(parsed);
    //var inheritancePattern = parsed.records[0].INFO.DNT;
    //var ret = {};
    //ret.parent0 = inheritancePattern.slice(0, 2);
    //ret.parent1 = inheritancePattern.slice(3, 5);
    //ret.child = inheritancePattern.slice(6, 8);
    //return ret;
  }

  module.exports.parseVCFText = parseVCFText;

}());
