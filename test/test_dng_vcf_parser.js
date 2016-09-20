var fs = require('fs');
var assert = require('assert');

var dngVCFParser = require('../src/dng_vcf_parser');

describe('dng vcf parser', function() {
  var vcfText = fs.readFileSync('dng.vcf', 'utf8');
  dngVCFParser.parseVCFText(vcfText);
});
