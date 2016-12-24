var fs = require('fs');
var assert = require('assert');

var dngVCFParser = require('../src/dng_vcf_parser');

describe('dng vcf parser', function() {
  it('read file', function() {
    var vcfText = fs.readFileSync('test/data/dng.vcf', 'utf8');
    dngVCFParser.parseVCFText(vcfText);
  });
});
