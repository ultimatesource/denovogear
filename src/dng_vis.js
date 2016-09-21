var vcfParser = require('./dng_vcf_parser');

jQuery.get('example_input_file.vcf', function(data) {
  d3.select('#text_box').text(data);
});

d3.select('#parse_button').on('click', parse);
d3.select('#vcf_file_input').on('change', function() {
  var selectedFile = document.getElementById('vcf_file_input').files[0];
  var reader = new FileReader();

  reader.onload = function(readerEvent) {
    d3.select('#text_box').text(reader.result);
  };
  reader.readAsText(selectedFile);
});

function parse() {
  // can't use d3 .text() method for this for some reason. Possibly d3 doesn't
  // work with textarea HTML elements
  var text = document.getElementById('text_box').value;
  var variantData = vcfParser.parseVCFText(text);
  console.log(variantData);
}
