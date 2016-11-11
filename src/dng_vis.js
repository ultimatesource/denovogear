var d3 = require('d3');
var jQuery = require('jquery');
//var graphml = require('graphml-js');
//var gl = require('graphlib');


jQuery.get('example_pedigree.ped', function(data) {
    d3.select('#text_box').text(data);
    main();
});

function main() {

  var width = 1280,
      height = 720;

  var svg = d3.select("#wrapper").append("svg")
      .attr("class", "mainSvg")
      .attr("width", width)
      .attr("height", height)
    .append("g");

  d3.select('#process_button').on('click', serverPedigreeAndLayout);
  d3.select('#pedigree_file_input').on('change', updateFile);

  serverPedigreeAndLayout();


  //var simulation = d3.forceSimulation()
  //  .force("link", d3.forceLink()
  //    .id(function(d) { return d.id; })
  //    .distance(function(d) { return 30; }))
  //  .force("charge", d3.forceManyBody()
  //    .strength(function(d) { return -50; }))
  //  .force("center", d3.forceCenter(width / 2, height / 2));


  //var parser = new graphml.GraphMLParser();

  //parser.parse(text, function (err, inputGraph) {
  //  var graph = convertToGraphlib(inputGraph);

  //  computeGenerations(graph);


  //  var graphmlGraph = convertFromGraphlib(graph);

  //  console.dir(graphmlGraph);
  //  doLayout(graphmlGraph, pedData);

  //  for (var node of graphmlGraph.nodes) {
  //    if (node.attributes.generation !== undefined) {
  //      node.fy = (node.attributes.generation * 200) - 180;
  //    }
  //  }

  //  var link = svg
  //    .append("g")
  //      .attr("class", "links")
  //    .selectAll("line")
  //    .data(graphmlGraph.edges)
  //    .enter().append("line")
  //      .attr("stroke-width", 1);


  //  //simulation
  //  //  .nodes(graphmlGraph.nodes)
  //  //  .on("tick", ticked);

  //  //simulation.force("link")
  //  //  .links(graphmlGraph.edges);

  //  function ticked() {
  //    link
  //      .attr("x1", function(d) { return d.source.x; })
  //      .attr("y1", function(d) { return d.source.y; })
  //      .attr("x2", function(d) { return d.target.x; })
  //      .attr("y2", function(d) { return d.target.y; });

  //    node.attr("transform", function(d) {
  //      return "translate(" + d.x + "," + d.y + ")";
  //    });
  //  }
  //});

  //function dragstarted(d) {
  //  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
  //  d.fx = d.x;

  //  if (!hasGeneration(d)) {
  //    d.fy = d.y;
  //  }
  //}

  //function dragged(d) {
  //  d.fx = d3.event.x;

  //  if (!hasGeneration(d)) {
  //    d.fy = d3.event.y;
  //  }
  //}

  //function dragended(d) {
  //  if (!d3.event.active) simulation.alphaTarget(0);
  //  d.fx = null;

  //  if (!hasGeneration(d)) {
  //    d.fy = null;
  //  }
  //}

  function serverPedigreeAndLayout() {
    var pedigreeText = document.getElementById('text_box').value;
    var pedigreeUploadData = { text: pedigreeText };
    jQuery.ajax('/pedigree_and_layout',
      { 
        type: 'POST',
        data: JSON.stringify(pedigreeUploadData),
        contentType: 'application/json',
        success: gotData
      });

    function gotData(jsonData) {
      var data = JSON.parse(jsonData);

      var ret = processPedigree(data);
      var nodes = ret.nodes;
      var links = ret.links;

      console.log(nodes);
      
      var link = svg
        .append("g")
          .attr("class", "links")
        .selectAll("line")
        .data(links)
        .enter().append("line")
          .attr("stroke-width", 1);
      link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });


      var node = svg
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
            if (d.type === 'male') {
              return d3.symbolSquare;
            }
            else if (d.type === 'female') {
              return d3.symbolCircle;
            }
            else {
              return d3.symbolTriangle;
            }
          })
          .size(500))
        .attr("fill", function(d) { 
          if (d.type === 'male') {
            return "SteelBlue";
          }
          else if (d.type === 'female') {
            return "Tomato";
          }
          else {
            return "black";
          }
        })
        .attr('opacity', function(d) {
          if (d.type === 'marriage') {
            return 0;
          }
          else {
            return 1;
          }
        });

      node.append("text")
        .attr("dx", 20)
        .attr("dy", ".35em")
        //.text(function(d) { return d.attributes.label });
        .text(function(d) { return d.dataNode.id });
      
      node.attr("transform", function(d) {
          return "translate(" + d.x + "," + d.y + ")";
      });

    }

  }

  function processPedigree(data) {
    var layout = data.layout;
    var nodes = [];
    var links = [];

    // build person nodes
    layout.nid.forEach(function(row, rowIdx) {
      for (var colIdx = 0; colIdx < layout.n[rowIdx]; colIdx++) {
        var id = row[colIdx];
        var node = {};
        node.dataNode = newNode(id);
        node.x = 100 + (100 * layout.pos[rowIdx][colIdx]);
        node.y = 100 * (rowIdx + 1);

        node.type = getGender(data.pedigree, id);

        // TODO: such a hack. remove
        node.rowIdx = rowIdx;
        node.colIdx = colIdx;

        nodes.push(node);
      }
    });

    // build marriage nodes
    nodes.forEach(function(node, index) {
      if (layout.spouse[node.rowIdx][node.colIdx] === 1) {
        var spouseNode = nodes[index + 1];
        var marriageNode = {};
        marriageNode.x = halfwayBetween(node, spouseNode);
        marriageNode.y = node.y;
        marriageNode.type = "marriage";
        marriageNode.dataNode = {};
        nodes.push(marriageNode);

        var marriageLink = {};
        marriageLink.source = node;
        marriageLink.target = spouseNode;
        links.push(marriageLink);

      }

      // TODO: such a hack. remove
      delete node.rowIdx;
      delete node.colIdx;
    });

    return { nodes: nodes, links: links };
  }

  function newNode(id) {
    return { id: id };
  }

  function halfwayBetween(nodeA, nodeB) {
    if (nodeB.x > nodeA.x) {
      return nodeA.x + ((nodeB.x - nodeA.x) / 2);
    }
    else {
      return nodeB.x + ((nodeA.x - nodeB.x) / 2);
    }
  }

  function getGender(pedigree, id) {
    for (var i in pedigree.id) {
      if (pedigree.id[i] === id) {
        return pedigree.sex[i];
      }
    }
  }

}

function updateFile() {
  var selectedFile = document.getElementById('pedigree_file_input').files[0];
  var reader = new FileReader();

  reader.onload = function(readerEvent) {
    d3.select('#text_box').text(reader.result);
  };
  reader.readAsText(selectedFile);
}

//function doLayout(graph, pedData) {
//  var numLevels = pedData.n.length;
//  for (var i = 0; i < numLevels; i++) {
//    var row = pedData.nid[i];
//    console.log("new row");
//    for (var j = 0; j < row.length; j++) {
//      var id = row[j];
//      if (id !== 0) {
//        var node = findNodeById(graph.nodes, id);
//        node.x = 100 + (50 * pedData.pos[i][j]);
//        node.y = 100 * (i + 1);
//        console.log(node.x);
//        console.log(node.y);
//      }
//    }
//  }
//}
//
//function findNodeById(nodes, id) {
//  var foundNode;
//  nodes.forEach(function(node, index) {
//    if (node.id === 'n' + id) {
//      foundNode = node;
//    }
//  });
//  return foundNode;
//}

//function parse() {
//  // can't use d3 .text() method for this for some reason. Possibly d3
//  // doesn't work with textarea HTML elements
//  var text = document.getElementById('text_box').value;
//  var parser = new graphml.GraphMLParser();
//
//  parser.parse(text, function (err, result) {
//    console.dir(result);
//  });
//}
//function convertToGraphlib(graph) {
//  var g = new gl.Graph({ directed: false });
//
//  for (node of graph.nodes) {
//    g.setNode(node.id, node.attributes);
//  }
//
//  for (edge of graph.edges) {
//    g.setEdge(edge.source, edge.target, edge.attributes);
//  }
//
//  return g;
//}

//function convertFromGraphlib(graph) {
//  var g = new graphml.Graph();
//
//  for (var nodeName of graph.nodes()) {
//    var newNode = new graphml.Node(nodeName);
//    newNode.attributes = graph.node(nodeName);
//    g.nodes.push(newNode);
//  }
//
//  for (var edgeObj of graph.edges()) {
//    var edgeName = edgeObj.v + "_" + edgeObj.w;
//    var newEdge = new graphml.Edge(edgeName, edgeObj.v, edgeObj.w);
//    newEdge.attributes = graph.edge(edgeObj);
//    g.edges.push(newEdge);
//  }
//
//  return g;
//}

//function computeGenerations(graph) {
//  var startNodeName = "n1";
//  var startNode = graph.node(startNodeName);
//  startNode.generation = 1;
//
//  removeEmptyNodes(graph);
//
//  computeGenerationsRecurse(graph, startNodeName);
//
//  normalizeGenerations(graph);
//}
//
//function computeGenerationsRecurse(graph, sourceName) {
//
//  var neighbors = graph.neighbors(sourceName);
//  var source = graph.node(sourceName);
//
//  for (var targetName of neighbors) {
//    var target = graph.node(targetName);
//
//    if (target.generation === undefined) {
//      var edge = graph.edge(sourceName, targetName);
//
//      if (edge.type === "Spousal") {
//        target.generation = source.generation;
//      }
//      else {
//        var sign;
//        if (isHigherGeneration(sourceName, targetName)) {
//          sign = 1;
//        }
//        else {
//          sign = -1;
//        }
//
//        if (edge.type === "Meiotic") {
//          target.generation = source.generation + (sign * 1);
//        }
//        else if (edge.type === "Mitotic") {
//          target.generation = source.generation + (sign * 0.2);
//        }
//        else if (edge.type === "Library") {
//          target.generation = source.generation + (sign * 0.2);
//        }
//        else {
//          throw new Error("Invalid Edge Type");
//        }
//      }
//
//      computeGenerationsRecurse(graph, targetName);
//    }
//  }
//}

//function isChildEdge(edge) {
//  return edge.attributes.type === 'Meiotic';
//}

//function isLibraryNode(node) {
//  if (node.attributes.label === undefined) {
//    return false;
//  }
//  return node.attributes.label.startsWith("LB");
//}

//function isSampleNode(node) {
//  if (node.attributes.label === undefined) {
//    return false;
//  }
//  return node.attributes.label.startsWith("SM");
//}

//function isSampleRootNode(node) {
//  return node.attributes.label === undefined;
//}

//function hasGeneration(node) {
//  return node.attributes.generation !== -1;
//}

//function isHigherGeneration(aName, bName) {
//  var aGen = Number(aName.substring(1));
//  var bGen = Number(bName.substring(1));
//  return bGen > aGen;
//}

//function normalizeGenerations(graph) {
//  var minGeneration = 1;
//
//  for (var nodeName of graph.nodes()) {
//    var node = graph.node(nodeName);
//
//    if (node.generation < minGeneration) {
//      // need to subtract 1 because the generations start at 1 istead of 0
//      minGeneration = node.generation - 1;
//    }
//  }
//
//  if (minGeneration < 1) {
//    for (var nodeName of graph.nodes()) {
//      var node = graph.node(nodeName);
//      node.generation = node.generation - minGeneration;
//    }
//  }
//}

//function removeEmptyNodes(graph) {
//  for (var nodeName of graph.nodes()) {
//    var node = graph.node(nodeName);
//    if (graph.neighbors(nodeName).length === 0) {
//      graph.removeNode(nodeName);
//    }
//  }
//}

//function numberOfGenerations(graph) {
//  var minGeneration = 0;
//  var maxGeneration = 0;
//
//  for (var nodeName of graph.nodes()) {
//    var node = graph.node(nodeName);
//
//    if (node.generation < minGeneration) {
//      // need to subtract 1 because the generations start at 1 istead of 0
//      minGeneration = node.generation;
//    }
//
//    if (node.generation > maxGeneration) {
//      // need to subtract 1 because the generations start at 1 istead of 0
//      maxGeneration = node.generation;
//    }
//  }
//
//  return Math.floor(maxGeneration - minGeneration);
//}

//function getId(node) {
//  return Number(node.id.substring(1));
//}

