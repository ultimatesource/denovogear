var d3 = require('d3');
var jQuery = require('jquery');

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
      console.log(jsonData);
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
        node.id = id;
        node.fatherId = data.pedigree.findex[oneToZeroBase(id)];
        node.motherId = data.pedigree.mindex[oneToZeroBase(id)];
        node.x = 410 + (80 * layout.pos[rowIdx][colIdx]);
        node.y = 100 + (100 * (rowIdx + 1));

        node.type = getGender(data.pedigree, id);

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
        var marriageNode = {};
        marriageNode.x = halfwayBetween(node, spouseNode);
        marriageNode.y = node.y;
        marriageNode.type = "marriage";
        marriageNode.dataNode = {};
        nodes.push(marriageNode);

        var marriageLeftLink = {};
        marriageLeftLink.type = 'spouse';
        marriageLeftLink.source = node;
        marriageLeftLink.target = marriageNode;
        links.push(marriageLeftLink);

        var marriageRightLink = {};
        marriageRightLink.type = 'spouse';
        marriageRightLink.source = spouseNode;
        marriageRightLink.target = marriageNode;
        links.push(marriageRightLink);

        var kids = getAllKids(data, nodes, node, spouseNode);

        for (var kid of kids) {
          var childLink = {};
          childLink.type = 'child';
          childLink.source = kid;
          childLink.target = marriageNode;
          links.push(childLink);
        }

      }

      // TODO: such a hack. remove
      delete node.rowIdx;
      delete node.colIdx;
    });

    return { nodes: nodes, links: links };
  }

  function oneToZeroBase(index) {
    return index - 1;
  }

  function newNode(id) {
    return { id: id };
  }

  function getAllKids(data, nodes, nodeA, nodeB) {
    var father;
    var mother;
    if (nodeA.type === 'male') {
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
        if (data.pedigree.findex[oneToZeroBase(node.id)] === father.id &&
            data.pedigree.mindex[oneToZeroBase(node.id)] === mother.id) {
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
