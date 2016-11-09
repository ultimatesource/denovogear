var express = require('express');
var bodyParser = require('body-parser');
var app = express();
var exec = require('child_process').exec;

app.use(express.static('dist'));
//app.use(bodyParser.urlencoded({ extended: false }));
app.use(bodyParser.json());

app.post('/pedigree_and_layout', function(req, res) {
  var r = exec('Rscript src/pedigree_and_layout.R', function(error, stdout) {
    res.end(stdout);
  });

  r.stdin.write(req.body.text);
  r.stdin.end();
});

app.listen(3000, function () {
    console.log('Hi there');
});
