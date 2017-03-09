// eslint exceptions
//
/* global vcf */
/* exported vcfParser */

var vcfParser = (function() {
  "use strict";

  function parseVCFText(text) {
    var parsed = vcf.parser()(text);
    parsed.header.contig[0].length = +parsed.header.contig[0].length;
    return parsed;
  }

  return {
    parseVCFText: parseVCFText
  };

}());
