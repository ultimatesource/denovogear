
main();

function main() {

  /*PEDIGREE_FILE_TEXT_PLACEHOLDER*/
  /*LAYOUT_DATA_PLACEHOLDER*/
  /*DNG_VCF_DATA_PLACEHOLDER*/

  var idText = d3.select('#id_display');
  var pedGraph = null;
  var activeNode = null;
  
  serverPedigreeAndLayout(function(nodes, links) {
    dngOverlay()
    visuals.doVisuals(nodes, links);
  });

  function dngOverlay() {
    var vcfData = vcfParser.parseVCFText(dngOutputFileText);

    var mutationLocation = vcfData.records[0].INFO.DNL;
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

    gotLayoutData();

    function gotLayoutData(jsonData) {

      var pedigreeData = pedParser.parsePedigreeFile(pedigreeFileText);

      var ret = processPedigree(layoutData, pedigreeData);
      var nodes = ret.nodes;
      var links = ret.links;

      callback(nodes, links);
    }
  }

  function processPedigree(data, pedigreeData) {

    pedGraph = buildGraphFromPedigree(pedigreeData);

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
    marriageNode.x = utils.halfwayBetween(spouseA.x, spouseB.x);
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
}
