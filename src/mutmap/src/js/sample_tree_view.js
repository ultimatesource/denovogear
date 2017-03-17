// eslint exceptions
//
/* global d3 */
/* exported sampleTreeView */

var sampleTreeView = (function(d3) {

  function createSampleTree() {

    var node;

    function my(selection) {

      var root = d3.hierarchy(node.dataNode.data.sampleIds);
      var treeWidth = 160;
      var treeHeight = 50;
      var cluster = d3.cluster().size([treeWidth, treeHeight]);
      cluster(root);

      var rootNode = root.descendants()[0];

      var sampleTree = selection.append("g")
          .attr("class", "sampleTree")
          .attr("transform", function() {
            // center the tree with the tree's root node overlapping the
            // current node
            return svgTranslateString(-rootNode.x, -rootNode.y);
          })
          .attr("visibility", "hidden");

      sampleTree.selectAll("sampleTreeLink")
          .data(root.descendants().slice(1))
        .enter().append("g")
          .attr("class", "sampleTreeLink")
        .append("line")
          .attr("stroke", "#999")
          .attr("x1", function(d) { return d.x; })
          .attr("y1", function(d) { return d.y; })
          .attr("x2", function(d) { return d.parent.x; })
          .attr("y2", function(d) { return d.parent.y; });

      var sampleTreeNodes = sampleTree.selectAll("sampleTreeNode")
          .data(root.descendants().slice(1))
        .enter().append("g")
          .attr("class", "sampleTreeNode")
          .attr("transform", function(d) {
            return svgTranslateString(d.x, d.y);
          });
        
      sampleTreeNodes.append("circle")
          .attr("r", 5)
          .attr("fill", "green");

      sampleTreeNodes.append("text")
          .attr("dx", 10)
          .attr("dy", ".35em")
          .text(function(d) {
            return d.data.name;
          });
    }

    my.node = function(value) {
      if (!arguments.length) return node;
      node = value;
      return my;
    };

    return my;
  }

  function svgTranslateString(x, y) {
    return "translate(" + x + "," + y + ")";
  }

  return {
    createSampleTree: createSampleTree
  };

}(d3));
