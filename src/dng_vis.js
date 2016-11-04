var d3 = require('d3');
var jQuery = require('jquery');
var graphml = require('graphml-js');
var gl = require('graphlib');

var pedData = {
  "n": [4, 5, 3, 1],
  "nid": [
    [1, 2, 4, 5, 0],
    [3, 6, 10, 9, 6],
    [7, 8, 11, 0, 0],
    [12, 0, 0, 0, 0]
  ],
  "pos": [
    [-2.4538e-16, 1, 4, 5, 0],
    [0.5, 1.5, 2.5, 3.5, 4.5],
    [0.7929, 1.7929, 2.7929, 0, 0],
    [2.2929, 0, 0, 0, 0]
  ],
  "fam": [
    [0, 0, 0, 0, 0],
    [1, 0, 0, 0, 3],
    [1, 1, 3, 0, 0],
    [2, 0, 0, 0, 0]
  ],
  "spouse": [
    [1, 0, 1, 0, 0],
    [1, 0, 1, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0]
  ]
};

jQuery.get('example_graph.graphml', function(data) {
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

  var simulation = d3.forceSimulation()
    .force("link", d3.forceLink()
      .id(function(d) { return d.id; })
      .distance(function(d) { return 30; }))
      //.strength(function(d) { return 3; }))
    .force("charge", d3.forceManyBody()
      .strength(function(d) { return -50; }))
    .force("center", d3.forceCenter(width / 2, height / 2));

    var text = document.getElementById('text_box').value;
    var parser = new graphml.GraphMLParser();

  parser.parse(text, function (err, inputGraph) {
    var graph = convertToGraphlib(inputGraph);

    computeGenerations(graph);


    var graphmlGraph = convertFromGraphlib(graph);

    console.dir(graphmlGraph);
    doLayout(graphmlGraph, pedData);

    for (var node of graphmlGraph.nodes) {
      if (node.attributes.generation !== undefined) {
        node.fy = (node.attributes.generation * 200) - 180;
      }
    }

    var link = svg
      .append("g")
        .attr("class", "links")
      .selectAll("line")
      .data(graphmlGraph.edges)
      .enter().append("line")
        .attr("stroke-width", 1);

    var node = svg
      .append("g")
        .attr("class", "nodes")
      .selectAll(".node")
      .data(graphmlGraph.nodes)
      .enter()
      .append("g")
        .attr("class", "node")
        .call(d3.drag()
          .on("start", dragstarted)
          .on("drag", dragged)
          .on("end", dragended));

    node.append("circle")
      .attr("r", 8)
      .attr("fill", function(d) { 
        if (isLibraryNode(d)) {
          return "Tomato";
        } 
        else if (isSampleNode(d)) {
          return "Grey";
        }
        else if (isSampleRootNode(d)) {
          return "DarkSeaGreen";
        }
        else {
          return "SteelBlue"; 
        }
      });

    node.append("text")
      .attr("dx", 12)
      .attr("dy", ".35em")
      //.text(function(d) { return d.attributes.label });
      .text(function(d) { return d.id });
    
    node.attr("transform", function(d) {
        return "translate(" + d.x + "," + d.y + ")";
    });

    //simulation
    //  .nodes(graphmlGraph.nodes)
    //  .on("tick", ticked);

    //simulation.force("link")
    //  .links(graphmlGraph.edges);

    function ticked() {
      link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

      node.attr("transform", function(d) {
        return "translate(" + d.x + "," + d.y + ")";
      });
    }
  });

  function dragstarted(d) {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;

    if (!hasGeneration(d)) {
      d.fy = d.y;
    }
  }

  function dragged(d) {
    d.fx = d3.event.x;

    if (!hasGeneration(d)) {
      d.fy = d3.event.y;
    }
  }

  function dragended(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = null;

    if (!hasGeneration(d)) {
      d.fy = null;
    }
  }
}

function parse() {
  // can't use d3 .text() method for this for some reason. Possibly d3
  // doesn't work with textarea HTML elements
  var text = document.getElementById('text_box').value;
  var parser = new graphml.GraphMLParser();

  parser.parse(text, function (err, result) {
    console.dir(result);
  });
}

d3.select('#parse_button').on('click', parse);
d3.select('#graph_file_input').on('change', updateGraphFile);

function updateGraphFile() {
  var selectedFile = document.getElementById('graph_file_input').files[0];
  var reader = new FileReader();

  reader.onload = function(readerEvent) {
    d3.select('#text_box').text(reader.result);
  };
  reader.readAsText(selectedFile);
}

function convertToGraphlib(graph) {
  var g = new gl.Graph({ directed: false });

  for (node of graph.nodes) {
    g.setNode(node.id, node.attributes);
  }

  for (edge of graph.edges) {
    g.setEdge(edge.source, edge.target, edge.attributes);
  }

  return g;
}

function convertFromGraphlib(graph) {
  var g = new graphml.Graph();

  for (var nodeName of graph.nodes()) {
    var newNode = new graphml.Node(nodeName);
    newNode.attributes = graph.node(nodeName);
    g.nodes.push(newNode);
  }

  for (var edgeObj of graph.edges()) {
    var edgeName = edgeObj.v + "_" + edgeObj.w;
    var newEdge = new graphml.Edge(edgeName, edgeObj.v, edgeObj.w);
    newEdge.attributes = graph.edge(edgeObj);
    g.edges.push(newEdge);
  }

  return g;
}

function computeGenerations(graph) {
  var startNodeName = "n1";
  var startNode = graph.node(startNodeName);
  startNode.generation = 1;

  removeEmptyNodes(graph);

  computeGenerationsRecurse(graph, startNodeName);

  normalizeGenerations(graph);
}

function computeGenerationsRecurse(graph, sourceName) {

  var neighbors = graph.neighbors(sourceName);
  var source = graph.node(sourceName);

  for (var targetName of neighbors) {
    var target = graph.node(targetName);

    if (target.generation === undefined) {
      var edge = graph.edge(sourceName, targetName);

      if (edge.type === "Spousal") {
        target.generation = source.generation;
      }
      else {
        var sign;
        if (isHigherGeneration(sourceName, targetName)) {
          sign = 1;
        }
        else {
          sign = -1;
        }

        if (edge.type === "Meiotic") {
          target.generation = source.generation + (sign * 1);
        }
        else if (edge.type === "Mitotic") {
          target.generation = source.generation + (sign * 0.2);
        }
        else if (edge.type === "Library") {
          target.generation = source.generation + (sign * 0.2);
        }
        else {
          throw new Error("Invalid Edge Type");
        }
      }

      computeGenerationsRecurse(graph, targetName);
    }
  }
}

function isChildEdge(edge) {
  return edge.attributes.type === 'Meiotic';
}

function isLibraryNode(node) {
  if (node.attributes.label === undefined) {
    return false;
  }
  return node.attributes.label.startsWith("LB");
}

function isSampleNode(node) {
  if (node.attributes.label === undefined) {
    return false;
  }
  return node.attributes.label.startsWith("SM");
}

function isSampleRootNode(node) {
  return node.attributes.label === undefined;
}

function hasGeneration(node) {
  return node.attributes.generation !== -1;
}

function isHigherGeneration(aName, bName) {
  var aGen = Number(aName.substring(1));
  var bGen = Number(bName.substring(1));
  return bGen > aGen;
}

function normalizeGenerations(graph) {
  var minGeneration = 1;

  for (var nodeName of graph.nodes()) {
    var node = graph.node(nodeName);

    if (node.generation < minGeneration) {
      // need to subtract 1 because the generations start at 1 istead of 0
      minGeneration = node.generation - 1;
    }
  }

  if (minGeneration < 1) {
    for (var nodeName of graph.nodes()) {
      var node = graph.node(nodeName);
      node.generation = node.generation - minGeneration;
    }
  }
}

function removeEmptyNodes(graph) {
  for (var nodeName of graph.nodes()) {
    var node = graph.node(nodeName);
    if (graph.neighbors(nodeName).length === 0) {
      graph.removeNode(nodeName);
    }
  }
}

function numberOfGenerations(graph) {
  var minGeneration = 0;
  var maxGeneration = 0;

  for (var nodeName of graph.nodes()) {
    var node = graph.node(nodeName);

    if (node.generation < minGeneration) {
      // need to subtract 1 because the generations start at 1 istead of 0
      minGeneration = node.generation;
    }

    if (node.generation > maxGeneration) {
      // need to subtract 1 because the generations start at 1 istead of 0
      maxGeneration = node.generation;
    }
  }

  return Math.floor(maxGeneration - minGeneration);
}

function doLayout(graph, pedData) {
  var numLevels = pedData.n.length;
  for (var i = 0; i < numLevels; i++) {
    var row = pedData.nid[i];
    console.log("new row");
    for (var j = 0; j < row.length; j++) {
      var id = row[j];
      if (id !== 0) {
        var node = findNodeById(graph.nodes, id);
        node.x = 100 + (50 * pedData.pos[i][j]);
        node.y = 100 * (i + 1);
        console.log(node.x);
        console.log(node.y);
      }
    }
  }
}

function getId(node) {
  return Number(node.id.substring(1));
}

function findNodeById(nodes, id) {
  var foundNode;
  nodes.forEach(function(node, index) {
    if (node.id === 'n' + id) {
      foundNode = node;
    }
  });
  return foundNode;
}
