
var vcfParserNew = (function() {
  "use strict";

  var optionsManager = utils.createOptionsManager();

  var FORMAT_INDEX = 8;
  var FIRST_SAMPLE_INDEX = 9;

  function VcfParser(options) {

  }

  VcfParser.prototype.parse = function(options) {

    optionsManager.checkOptions({
      requiredOptions: ['vcfText'],
      providedOptions: options
    });

    var vcfText = options.vcfText;

    this._vcfData = {
      header: {
        contig: [],
        sampleNames: [],
      },
      records: [],
    };

    var lines = vcfText.split("\n");
    lines.forEach(function(line, index) {
      if (line.startsWith("##contig")) {
        this._parseContig(line);
      }
      else if (line.startsWith("#CHROM")) {
        this._parseHeader(line);
      }
      else if (!line.startsWith("#") && line.length > 0) {

        var columns = this._parseTSVLine(line);
        var chrom = this._parseChromosome(columns);
        
        var parser = this;
        var cache = {};

        var record = {
          // TODO: Cache the CHROM just like the other properties. For some
          // reason this is currently very slow. I'm concerned there may be a
          // more subtle bug occurring.
          CHROM: chrom,
          get POS() {
            if (cache.POS === undefined) {
              cache.POS = parser._parsePosition(columns);
            }
            return cache.POS;
          },
          get INFO() {
            if (cache.INFO === undefined) {
              cache.INFO = parser._parseInfoColumn(columns);
            }
            return cache.INFO;
          },
        };

        this._vcfData.header.sampleNames.forEach(function(sampleName, i) {
          var sampleIndex = i + FIRST_SAMPLE_INDEX;
          Object.defineProperty(record, sampleName, {
            get: function() {
              if (cache[sampleName] === undefined) {
                var format = parser._parseFormat(columns[FORMAT_INDEX]);
                cache[sampleName] =
                  parser._parseSampleInfo(columns[sampleIndex], format);
              }
              return cache[sampleName];
            }
          });
        }, this);

        this._vcfData.records.push(record);
      }
    }, this);

    return this._vcfData;
  }

  VcfParser.create = function(options) {
    return new VcfParser(options);
  };

  VcfParser.prototype._parseContig = function(line) {
    var start = line.indexOf("<");
    var comma = line.indexOf(",");
    var end = line.indexOf(">");

    var idPair = this._parsePair(line.slice(start + 1, comma));
    var id = idPair.value;
    var lengthPair = this._parsePair(line.slice(comma, end));
    //console.log(lengthPair);
    var length = Number(lengthPair.value);

    this._vcfData.header.contig.push({ ID: id, length: length });
  };

  VcfParser.prototype._parseHeader = function(line) {
    var columns = this._parseTSVLine(line);
    for (var i = FIRST_SAMPLE_INDEX; i < columns.length; i++) {
      this._vcfData.header.sampleNames.push(columns[i]);
    }
  };

  VcfParser.prototype._parseTSVLine = function(line) {
    var columns = line.split("\t");
    return columns;
  };

  VcfParser.prototype._parseInfoColumn = function(columns) {

    var infoColumn = columns[7];
    var pairs = infoColumn.split(';');

    var info = {};
    pairs.forEach(function(pairString) {
      var pair = this._parsePair(pairString);
      info[pair.key] = pair.value;
    }, this);

    return info;
  };

  VcfParser.prototype._parseChromosome = function(columns) {
    return columns[0];
  };

  VcfParser.prototype._parsePosition = function(columns) {
    return Number(columns[1]);
  };

  VcfParser.prototype._parseSamples = function(columns) {

    var samples = {};
    var format = this._parseFormat(columns[FORMAT_INDEX]);
    for (var i = FIRST_SAMPLE_INDEX; i < columns.length; i++) {
      var sampleName =
        this._vcfData.header.sampleNames[i - FIRST_SAMPLE_INDEX];
      samples[sampleName] = this._parseSampleInfo(columns[i], format);
    }

    return samples;

  };

  VcfParser.prototype._parseSample = function(column, format) {
      var sampleName =
        this._vcfData.header.sampleNames[i - FIRST_SAMPLE_INDEX];
      samples[sampleName] = this._parseSampleInfo(columns, format);
  };

  VcfParser.prototype._parseSampleInfo = function(column, format) {

    var items = column.split(":");
    var sampleInfo = {};
    for (var i = 0; i < format.length; i++) {

      if (format[i] === "GT") {
        sampleInfo[format[i]] = items[i];
      }
      else if (items[i].includes(",")) {
        sampleInfo[format[i]] = items[i].split(",");

        for (var j = 0; j < sampleInfo[format[i]].length; j++ ) {
          sampleInfo[format[i]][j] = Number(sampleInfo[format[i]][j]);
        }
      }
      else {
        sampleInfo[format[i]] = Number(items[i]);
      }
    }

    return sampleInfo;
  };

  VcfParser.prototype._parseFormat = function(column) {
    return column.split(":");
  };

  VcfParser.prototype._parsePair = function(pair) {
    var pairArray = pair.split('=');
    return {
      key: pairArray[0],
      value: pairArray[1]
    };
  }

  return {
    VcfParser: VcfParser
  };


}());
