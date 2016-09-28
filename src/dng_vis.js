var parseString = require('xml2js').parseString;

jQuery.get('example_graph.graphml', function(data) {
  d3.select('#text_box').text(data);
});

d3.select('#parse_button').on('click', parse);
d3.select('#graph_file_input').on('change', function() {
  var selectedFile = document.getElementById('graph_file_input').files[0];
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
  parseString(text, function (err, result) {
      console.dir(result);
  });
}
