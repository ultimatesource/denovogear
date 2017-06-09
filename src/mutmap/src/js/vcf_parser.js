
var vcfParserNew = (function() {
  "use strict";

  var optionsManager = utils.createOptionsManager();

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
      },
      records: [],
    };

    var lines = vcfText.split("\n");
    lines.forEach(function(line, index) {
      if (line.startsWith("##contig")) {
        this._parseContig(line);
      }
      else if (!line.startsWith("#") && line.length > 0) {


        var columns = this._parseDataLine(line);
        var chrom = this._parseChromosome(columns);
        var pos = this._parsePosition(columns);
        //var info = this._parseInfoColumn(data[7]);
        //var dnl = parseDNL(info);
        
        var record = {
          chrom: chrom,
          pos: pos,
        };

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

    this._vcfData.header.contig.push({ id: id, length: length });
  };

  VcfParser.prototype._parseDataLine = function(line) {
    var columns = line.split("\t");
    return columns;
  };

  VcfParser.prototype._parseInfoColumn = function(column) {
    return column.split(';');
  };

  VcfParser.prototype._parseChromosome = function(columns) {
    return columns[0];
  };

  VcfParser.prototype._parsePosition = function(columns) {
    return Number(columns[1]);
  };

  VcfParser.prototype._parseDNL = function(info) {
    return this._parsePair(info[6]);
  }

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
