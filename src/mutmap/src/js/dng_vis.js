// eslint exceptions
//
/* global pedParser */
/* global visuals */
/* global vcfParser */
/* global pedigr */
/* global utils */

/* global pedigreeFileText */
/* global layoutData */
/* global dngOutputFileText */

(function() {
  "use strict";

  // The placeholder tags below will be replaced with the correct data objects
  // by the build system. See mutmap/tools/layout_and_template.R
  //
  // provides pedigreeFileText
  /*PEDIGREE_FILE_TEXT_PLACEHOLDER*/
  // provides layoutData
  /*LAYOUT_DATA_PLACEHOLDER*/
  // provides dngOutputFileText
  /*DNG_VCF_DATA_PLACEHOLDER*/

  
  var pedigreeData = pedParser.parsePedigreeFile(pedigreeFileText);

  var pedGraph = buildGraphFromPedigree(pedigreeData);

  console.log(pedGraph);

  var kinshipPedigreeData = layoutData;

  var graphData = processPedigree(kinshipPedigreeData);

  dngOverlay();
  visuals.doVisuals(graphData);

  window.addEventListener("resize", function() {
    visuals.doVisuals(graphData);
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

      vcfData.header.sampleNames.forEach(function(sampleName) {
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
      });
    }
    else {
      alert("No mutation found!");
    }
  }

  function processPedigree(kinshipPedigreeData) {

    var layout = kinshipPedigreeData.layout;
    var nodes = [];
    var links = [];

    // build person nodes
    layout.nid.forEach(function(row, rowIdx) {
      for (var colIdx = 0; colIdx < layout.n[rowIdx]; colIdx++) {

        var id = row[colIdx];

        var node = {};
        node.type = "person";
        node.dataNode = pedGraph.getPerson(id);
        node.x = 120 * layout.pos[rowIdx][colIdx];
        node.y = 120 * rowIdx;

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

        var children = getAllChildren(kinshipPedigreeData, nodes, node,
                                      spouseNode);

        var encountered = [];

        children.forEach(function(childNode) {
          var index = encountered.findIndex(function(element) {
            return element.dataNode === childNode.dataNode;
          });

          if (index === -1) {
            // TODO: this is duplicated below. fix it
            var childLink = createChildLink(childNode, marriageNode);
            var parentageLink = marriage.addChild(childNode.dataNode);
            childLink.dataLink = parentageLink;
            links.push(childLink);
            encountered.push(childNode);
          }
          else {
            var oldOne = encountered[index];
            var distanceOldOneToParents = distanceBetweenNodes(oldOne,
              marriageNode);
            var distanceCurrentToParents = distanceBetweenNodes(childNode,
              marriageNode);
            
            if (distanceOldOneToParents > distanceCurrentToParents) {
              var index = links.findIndex(function(element) {
                return element.type === "child" &&
                  element.source.dataNode === childNode.dataNode;
              });
              links.splice(index, 1);

              var childLink = createChildLink(childNode, marriageNode);
              var parentageLink = marriage.addChild(childNode.dataNode);
              childLink.dataLink = parentageLink;
              links.push(childLink);
            }

            var duplicateLink = createDuplicateLink(childNode, oldOne);
            links.push(duplicateLink);
          }
        });

      }

      // TODO: such a hack. remove
      delete node.rowIdx;
      delete node.colIdx;
    });

    return { nodes: nodes, links: links };
  }

  function buildGraphFromPedigree(pedigreeData) {
    var pedGraph = pedigr.PedigreeGraph.createGraph();

    pedigreeData.forEach(function(individual) {
      var person = pedigr.PersonBuilder
        .createPersonBuilder(individual.individualId)
          .sex(individual.sex)
          .data({ sampleIds: individual.sampleIds })
          .build();
      pedGraph.addPerson(person);
    });

    return pedGraph;
  }

  // TODO this and findMatchingSampleNode have almost the same logic. Find a
  // way to extract the duplication
  function findOwnerNode(sampleName) {
    var strippedName = getStrippedName(sampleName);
    var persons = pedGraph.getPersons();
    for (var index = 0; index < persons.length; index++) {
      var person = persons[index];
      var sampleNode = findInTree(person.data.sampleIds, strippedName);
      if (sampleNode !== undefined) {
        return person;
      }
    }
    return undefined;
  }

  function findMatchingSampleNode(sampleName) {
    var strippedName = getStrippedName(sampleName);
    var persons = pedGraph.getPersons();
    for (var index = 0; index < persons.length; index++) {
      var person = persons[index];
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
        for (var index = 0; index < tree.children.length; index++) {
          var child = tree.children[index];
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
    var stripped = sampleName.slice(3, sampleName.indexOf(":"));
    return stripped;
  }

  function isPersonNode(sampleName) {
    return sampleName.startsWith("GL-");
  }

  function getIdFromSampleName(sampleName) {
    return sampleName.slice(3);
  }

  function oneToZeroBase(index) {
    return index - 1;
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
    marriageLink.type = "spouse";
    marriageLink.source = spouseNode;
    marriageLink.target = marriageNode;
    return marriageLink;
  }

  function createChildLink(childNode, marriageNode) {
    var childLink = {};
    childLink.type = "child";
    childLink.source = childNode;
    childLink.target = marriageNode;
    return childLink;
  }

  function createDuplicateLink(nodeA, nodeB) {
    var dupLink = {};
    dupLink.type = "duplicate";
    dupLink.source = nodeA;
    dupLink.target = nodeB;
    return dupLink;
  }

  function getAllChildren(kinshipPedigreeData, nodes, nodeA, nodeB) {
    var father;
    var mother;
    if (nodeA.dataNode.sex === "male") {
      father = nodeA;
      mother = nodeB;
    }
    else {
      father = nodeB;
      mother = nodeA;
    }

    var children = [];
    nodes.forEach(function(node) {
      if (node.type != "marriage") {
        if (kinshipPedigreeData.pedigree.findex[oneToZeroBase(node.dataNode.id)] ===
              father.dataNode.id &&
            kinshipPedigreeData.pedigree.mindex[oneToZeroBase(node.dataNode.id)] ===
              mother.dataNode.id) {
          children.push(node);
        }
      }
    });

    return children;
  }

  function distanceBetweenNodes(nodeA, nodeB) {
    return utils.distanceBetweenPoints(nodeA.x, nodeA.y, nodeB.x, nodeB.y);
  }
  
}());
