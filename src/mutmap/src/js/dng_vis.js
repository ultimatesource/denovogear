//var d3 = require('d3');
//var vcfParser = require('./dng_vcf_parser');
//var pedParser = require('./ped_parser');
//var pedigr = require('./pedigr');

//var pedigreeFileText;
//var dngOutputFileText;

//TODO: clean this up
//jQuery.get('example_pedigree.ped', function(pedData) {
//    pedigreeFileText = pedData;
//
//    jQuery.get('example_output.vcf', function(outputData) {
//      dngOutputFileText = outputData;
//      main();
//    });
//});

main();

function main() {

  /*PEDIGREE_FILE_TEXT_PLACEHOLDER*/
  /*LAYOUT_DATA_PLACEHOLDER*/
  /*DNG_VCF_DATA_PLACEHOLDER*/

  //dngOutputData = dngOutputData.replace('\\\n', function() { return '\n'; });
  //dngOutputData = dngOutputData.replace('\"', function() { return '"'; });

  //d3.select('#pedigree_file_input').on('change', updatePedigreeFile);
  //d3.select('#dng_output_file_input').on('change', updateDNGOutputFile);

  var idText = d3.select('#id_display');
  var pedGraph = null;
  var activeNode = null;
  
  serverPedigreeAndLayout(function(nodes, links) {
    dngOverlay()
    doVisuals(nodes, links);
  });

  function dngOverlay() {
    var vcfData = vcfParser.parseVCFText(dngOutputFileText);
    //console.log(vcfData);

    var mutationLocation = vcfData.records[0].INFO.DNL;
    console.log(mutationLocation);
    var owner = findOwnerNode(mutationLocation);

    if (owner !== undefined) {
      var ownerParentageLink = owner.getParentageLink();
      var parentageLinkData = {
        mutation: vcfData.records[0].INFO.DNT
      };
      ownerParentageLink.setData(parentageLinkData);

      for (var sampleName of vcfData.header.sampleNames) {
        var format = vcfData.records[0][sampleName];

        if (isPersonNode(sampleName)) {
          var id = getIdFromSampleName(sampleName);
          var personNode = pedGraph.getPerson(id);
          personNode.data.dngOutputData = format;
        }
        else {
          var sampleNode = findMatchingSampleNode(sampleName);
          sampleNode.dngOutputData = format;
        }
      }
    }
    else {
      alert("No mutation found!");
    }
  }

  function serverPedigreeAndLayout(callback) {

    //var pedigreeUploadData = { text: pedigreeFileText };
    //jQuery.ajax('/pedigree_and_layout',
    //  { 
    //    type: 'POST',
    //    data: JSON.stringify(pedigreeUploadData),
    //    contentType: 'application/json',
    //    success: gotLayoutData
    //  });
    gotLayoutData();

    function gotLayoutData(jsonData) {
      //console.log(jsonData);
      //var layoutData = JSON.parse(jsonData);

      var pedigreeData = pedParser.parsePedigreeFile(pedigreeFileText);

      var ret = processPedigree(layoutData, pedigreeData);
      var nodes = ret.nodes;
      var links = ret.links;

      //doVisuals(nodes, links);

      callback(nodes, links);
    }
  }

  function doVisuals(nodes, links) {
      var height = 500;

      d3.select("svg").remove();

      var zoom = d3.zoom()
        .on("zoom", zoomed);


      var chartWrapper= d3.select("#chart_wrapper")
      var dim = chartWrapper.node().getBoundingClientRect();

      var svg = chartWrapper.append("svg")
          .attr("width", dim.width)
          .attr("height", height);

      var container = svg.call(zoom)
        .append("g");

        
      var link = container
        .append("g")
          .attr("class", "links")
        .selectAll("line")
        .data(links)
        .enter()
        .append("g")
          .attr("class", "link");

      link.append("line")
        .attr("stroke-width", function(d) {
          if (linkHasData(d)) {
              return 5;
          }
          return 1;
        })
        .attr("stroke", function(d) {
          if (linkHasData(d)) {
            return "green";
          }

          return "#999";
        })
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

      link.append("text")
        .attr("dx", function(d) {
          return halfwayBetweenGeneric(d.source.x, d.target.x) - 45;
        })
        .attr("dy", function(d) {
          return halfwayBetweenGeneric(d.source.y, d.target.y);
        })
        .text(function(d) {
          if (linkHasData(d)) {
            return d.dataLink.data.mutation;
          }
        });

      function linkHasData(d) {
          return d.dataLink !== undefined && d.dataLink.data !== undefined;
      }


      var node = container
        .append("g")
          .attr("class", "nodes")
        .selectAll(".node")
        .data(nodes)
        .enter()
        .append("g")
          .attr("class", "node");

      node.append("path")
        .attr("d", d3.symbol()
          .type(function(d) {
            if (d.type === 'person') {
              if (d.dataNode.sex === 'male') {
                return d3.symbolSquare;
              }
              else if (d.dataNode.sex === 'female') {
                return d3.symbolCircle;
              }
            }
            else {
              return d3.symbolTriangle;
            }
          })
          .size(500))
        .attr("fill", fillColor)
        .attr('opacity', function(d) {
          if (d.type === 'marriage') {
            return 0;
          }
          else {
            return 1;
          }
        })
        .on("click", nodeClicked);

      node.append("text")
        .attr("dx", 20)
        .attr("dy", ".35em")
        .text(function(d) { 
          if (d.type !== 'marriage') {
            return d.dataNode.id
          }
        });
      
      node.attr("transform", function(d) {
          return "translate(" + d.x + "," + d.y + ")";
      });

      function zoomed() {
        container.attr("transform", d3.event.transform);
      }

      function nodeClicked(d) {
        if (d.type !== 'marriage') {
          d3.select(activeNode).style('fill', fillColor);
          activeNode = this;
          d3.select(this).style('fill', 'DarkSeaGreen');

          if (d.dataNode.data.dngOutputData !== undefined) {
            document.getElementById('id_display').value =
              d.dataNode.id;
          }
          else {
            document.getElementById('id_display').value = "";
          }

        }
      }

      function fillColor(d) {
        if (d.type === 'person') {
          if (d.dataNode.sex === 'male') {
            return "SteelBlue";
          }
          else if (d.dataNode.sex === 'female') {
            return "Tomato";
          }
        }
        else {
          return "black";
        }
      }
  }

  function processPedigree(data, pedigreeData) {
    //console.log(pedigreeData);

    pedGraph = buildGraphFromPedigree(pedigreeData);

    //console.log(pedGraph);

    var layout = data.layout;
    var nodes = [];
    var links = [];

    // build person nodes
    layout.nid.forEach(function(row, rowIdx) {
      for (var colIdx = 0; colIdx < layout.n[rowIdx]; colIdx++) {

        var id = row[colIdx];

        var node = {};
        node.type = 'person';
        node.dataNode = pedGraph.getPerson(id);
        node.x = 80 * layout.pos[rowIdx][colIdx];
        node.y = 100 * rowIdx;

        // TODO: such a hack. remove
        node.rowIdx = rowIdx;
        node.colIdx = colIdx;

        nodes.push(node);
      }
    });

    // build marriage nodes and links
    nodes.forEach(function(node, index) {
      if (layout.spouse[node.rowIdx][node.colIdx] === 1) {
        var spouseNode = nodes[index + 1];

        var marriageNode = createMarriageNode(node, spouseNode);
        nodes.push(marriageNode);

        links.push(createMarriageLink(node, marriageNode));
        links.push(createMarriageLink(spouseNode, marriageNode));

        var marriage = pedigr.MarriageBuilder.createMarriageBuilder()
          .spouse(node.dataNode)
          .spouse(spouseNode.dataNode)
          .build();

        pedGraph.addMarriage(marriage);

        var children = getAllKids(data, nodes, node, spouseNode);

        for (var childNode of children) {
          var childLink = createChildLink(childNode, marriageNode);
          var parentageLink = marriage.addChild(childNode.dataNode);
          childLink.dataLink = parentageLink;
          links.push(childLink);
        }

      }

      // TODO: such a hack. remove
      delete node.rowIdx;
      delete node.colIdx;
    });

    return { nodes: nodes, links: links };
  }

  function buildGraphFromPedigree(pedigreeData) {
    var pedGraph = pedigr.PedigreeGraph.createGraph();

    pedigreeData.forEach(function(person) {
      var person = pedigr.PersonBuilder
        .createPersonBuilder(person.individualId)
          .sex(person.sex)
          .data({ sampleIds: person.sampleIds })
          .build();
      pedGraph.addPerson(person);
    });

    return pedGraph;
  }

  function findOwnerNode(sampleName) {
    var strippedName = getStrippedName(sampleName);
    for (var person of pedGraph.getPersons()) {
      console.log(person.data.sampleIds);
      var sampleNode = findInTree(person.data.sampleIds, strippedName);
      if (sampleNode !== undefined) {
        return person;
      }
    }
    return undefined;
  }

  function findMatchingSampleNode(sampleName) {
    var strippedName = getStrippedName(sampleName);
    for (var person of pedGraph.getPersons()) {
      var sampleNode = findInTree(person.data.sampleIds, strippedName);
      if (sampleNode !== undefined) {
        return sampleNode;
      }
    }
    return undefined;
  }

  function findInTree(tree, sampleName) {

    if (tree.name === sampleName) {
      return tree;
    }

    if (tree.children !== undefined) {
      if (tree.children.length === 0) {
        return undefined;
      }
      else {
        for (child of tree.children) {
          var inChild = findInTree(child, sampleName);
          if (inChild !== undefined) {
            return inChild;
          }
        }
      }
    }

    return undefined;
  }

  function getStrippedName(sampleName) {
    var stripped = sampleName.slice(3, sampleName.indexOf(':'));
    return stripped;
  }

  function isPersonNode(sampleName) {
    return sampleName.startsWith('GL-');
  }

  function getIdFromSampleName(sampleName) {
    return sampleName.slice(3);
  }

  function oneToZeroBase(index) {
    return index - 1;
  }

  function newNode(id) {
    return { id: id };
  }

  function createMarriageNode(spouseA, spouseB) {
    var marriageNode = {};
    marriageNode.x = halfwayBetween(spouseA, spouseB);
    marriageNode.y = spouseA.y;
    marriageNode.type = "marriage";
    marriageNode.dataNode = {};
    return marriageNode;
  }

  function createMarriageLink(spouseNode, marriageNode) {
    var marriageLink = {};
    marriageLink.type = 'spouse';
    marriageLink.source = spouseNode;
    marriageLink.target = marriageNode;
    return marriageLink;
  }

  function createChildLink(childNode, marriageNode) {
    var childLink = {};
    childLink.type = 'child';
    childLink.source = childNode;
    childLink.target = marriageNode;
    return childLink;
  }

  function getAllKids(data, nodes, nodeA, nodeB) {
    var father;
    var mother;
    if (nodeA.dataNode.sex === 'male') {
      father = nodeA;
      mother = nodeB;
    }
    else {
      father = nodeB;
      mother = nodeA;
    }

    var kids = [];
    for (var node of nodes) {
      if (node.type != 'marriage') {
        if (data.pedigree.findex[oneToZeroBase(node.dataNode.id)] ===
              father.dataNode.id &&
            data.pedigree.mindex[oneToZeroBase(node.dataNode.id)] ===
              mother.dataNode.id) {
          kids.push(node);
        }
      }
    }

    return kids;
  }

  function halfwayBetween(nodeA, nodeB) {
    if (nodeB.x > nodeA.x) {
      return nodeA.x + ((nodeB.x - nodeA.x) / 2);
    }
    else {
      return nodeB.x + ((nodeA.x - nodeB.x) / 2);
    }
  }

  function halfwayBetweenGeneric(a, b) {
    if (a > b) {
      return a + ((b - a) / 2);
    }
    else {
      return b + ((a - b) / 2);
    }
  }

  //function updatePedigreeFile() {
  //  updateFile('pedigree_file_input', function(fileData) {
  //    pedigreeFileText = fileData;
  //    serverPedigreeAndLayout(function(nodes, links) {
  //      doVisuals(nodes, links);
  //    });
  //  });
  //}

  //function updateDNGOutputFile() {
  //  updateFile('dng_output_file_input', function(fileData) {
  //    serverPedigreeAndLayout(function(nodes, links) {
  //      dngOutputFileText = fileData;
  //      dngOverlay();
  //      doVisuals(nodes, links);
  //    });
  //  });
  //}

  //function updateFile(fileInputElementId, callback) {
  //  var selectedFile = document.getElementById(fileInputElementId).files[0];

  //  if (selectedFile !== undefined) {
  //    var reader = new FileReader();

  //    reader.onload = function(readerEvent) {
  //      callback(reader.result);
  //    };
  //    reader.readAsText(selectedFile);
  //  }
  //}
}
