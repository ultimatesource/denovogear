var fs = require('fs');
var vcf = require('vcf.js');

var contents = fs.readFileSync('dng.vcf', 'utf8');
console.log(contents);

var parsed = vcf.parser()(contents);

var inheritancePattern = parsed['records'][0]['INFO']['DNT'];
var parent0 = inheritancePattern.slice(0, 2);
var parent1 = inheritancePattern.slice(3, 5);
var child = inheritancePattern.slice(6, 8);
console.log(parent0);
console.log(parent1);
console.log(child);
//console.log(parsed['header']['FORMAT']);
