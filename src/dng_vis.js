var d3 = require('d3');
var jQuery = require('jquery');
var graphml = require('graphml-js');
var gl = require('graphlib');

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
    console.dir(inputGraph);

    var graph = convertToGraphlib(inputGraph);

    process(graph);
    var graphmlGraph = convertFromGraphlib(graph);
    console.dir(graphmlGraph);

    for (var node of graphmlGraph.nodes) {
      if (node.attributes.generation !== -1) {
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
          return "tomato";
        } 
        else if (isSampleNode(d)) {
          return "grey";
        }
        else {
          return "steelblue"; 
        }
      });

    node.append("text")
      .attr("dx", 12)
      .attr("dy", ".35em")
      //.text(function(d) { return d.attributes.label });
      .text(function(d) { return d.id });

    simulation
      .nodes(graphmlGraph.nodes)
      .on("tick", ticked);

    simulation.force("link")
      .links(graphmlGraph.edges);

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

function process(graph) {

  for (var sourceName of graph.nodes()) {
    var source = graph.node(sourceName);

    if (source.generation === undefined) {
      source.generation = 1;
    }

    var neighbors = graph.neighbors(sourceName);

    // remove any nodes with no edges
    if (neighbors.length === 0) {
      graph.removeNode(sourceName);
    }
    else {
      for (var targetName of neighbors) {
        var target = graph.node(targetName);

        if (target.generation === undefined) {
          var edge = graph.edge(sourceName, targetName);

          if (edge.type === "Meiotic") {
            target.generation = source.generation + 1;
          }
          else if (edge.type === "Spousal") {
            target.generation = source.generation;
          }
          else if (edge.type === "Mitotic") {
            target.generation = -1;
          }
          else if (edge.type === "Library") {
            target.generation = -1;
          }
          else {
            throw new Error("Fail");
          }
        }
      }
    }
  }
}

function isLibraryNode(node) {
  return node.attributes.label.startsWith("LB");
}

function isSampleNode(node) {
  return node.attributes.label.startsWith("SM");
}

function hasGeneration(node) {
  return node.attributes.generation !== -1;
}
