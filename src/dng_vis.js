var vcfParser = require('./dng_vcf_parser');

d3.select('#parse_button').on('click', parse);

function parse() {
  var text = document.getElementById('text_box').value;
  var variantData = vcfParser.parseVCFText(text);
  console.log(variantData);
}
