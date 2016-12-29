var pedParser = (function() {
  "use strict";

  function parsePedigreeFile(fileText) {

    var newickParser = require('biojs-io-newick');

    var entries = [];

    var lines = fileText.split('\n');
    for (var line of lines) {

      if (line.length === 0) {
        continue;
      }

      var row = line.split("\t");
      var entry = {};
      entry.familyId = Number(row[0]);
      entry.individualId = Number(row[1]);

      var fatherId = Number(row[2]);
      if (fatherId === 0) {
        entry.fatherId = undefined;
      }
      else {
        entry.fatherId = fatherId;
      }

      var motherId = Number(row[3]);
      if (motherId === 0) {
        entry.motherId = undefined;
      }
      else {
        entry.motherId = motherId;
      }

      var sex = Number(row[4]);
      if (sex === 1) {
        entry.sex = 'male';
      }
      else if (sex === 2) {
        entry.sex = 'female';
      }
      else {
        entry.sex = undefined;
      }

      var sampleIds = row[5];

      if (sampleIds !== undefined) {
        entry.sampleIds =
          newickParser.parse_newick(convertToValidNewick(sampleIds));
      }
      else {
        entry.sampleIds = {
          children: [],
          name: String(entry.individualId)
        };
      }

      entries.push(entry);
    }

    return entries;
  }

  function convertToValidNewick(row) {
    return "(" + row + ")";
  }

  return {
    parsePedigreeFile: parsePedigreeFile
  }

}());
