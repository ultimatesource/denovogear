var fs = require('fs');

var templateData = fs.readFileSync('index.html', { encoding: 'utf-8' });
var bundleData = fs.readFileSync('build/bundle.js', { encoding: 'utf-8' });

// Note: need to use function replacement argument because if you use a string,
// and the string contains '$', it has special replacement meaning.
var outData = templateData.replace('<!--BUNDLE_JS_PLACEHOLDER-->', function() {
  return bundleData;
});

fs.writeFileSync('dist/index.html', outData);
