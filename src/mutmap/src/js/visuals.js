var visuals = (function() {
  "use strict";
  
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
        return utils.halfwayBetween(d.source.x, d.target.x) - 45;
      })
      .attr("dy", function(d) {
        return utils.halfwayBetween(d.source.y, d.target.y);
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

  return { doVisuals: doVisuals };

}());
