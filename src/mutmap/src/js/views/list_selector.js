mutmap.ListSelectorView = (function(d3, utils) {
  "use strict";
  
  var ListSelectorView = Backbone.View.extend({
    
    initialize: function(options) {
      this.render(options);
    },

    render: function(options) {

      var model = this.model;

      var list = options.list;
      var selector = options.selector;
      var itemName = options.itemName;

      var currentChromIndex = 0;

      this._renderInto = this.el;

      var container = this._renderInto.append("div")
          .attr("class", "list-selector");

      container.append("div")
        .append("text")
          .text(currentItem);

      var prevButton = container.append("button")
          .attr("class", "btn btn-default")
          .on("click", function() {
            currentChromIndex--;
            if (currentChromIndex < 0) {
              currentChromIndex = list.length - 1;
            }

            model.set('index', currentChromIndex);
          });
      prevButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-left");
      prevButton.append("span")
          .text(" Previous");

      var dropdown = container.append("div")
          .attr("class", "dropdown")
      dropdown.append("button")
          .attr("class", "btn btn-default dropdown-toggle")
          .attr("id", "dropdownMenu1")
          .attr("type", "button")
          .attr("data-toggle", "dropdown")
          .text("Select " + itemName);

      dropdown.append("ul")
          .attr("class", "dropdown-menu")
        .selectAll(".mutation")
          .data(list)
        .enter().append("li")
          .attr("class", "mutation")
        .append("a")
          .attr("href", "#")
          .text(function(d) {
            if (selector !== undefined) {
              return d[selector];
            }
            else {
              return d;
            }
          })
          .on("click", function(d, i) {
            currentChromIndex = i;
            model.set('index', currentChromIndex);
          });

        var nextButton = container.append("button")
            .attr("class", "btn btn-default")
            .on("click", function() {
              currentChromIndex++;
              if (currentChromIndex === list.length) {
                currentChromIndex = 0;
              }

              model.set('index', currentChromIndex);
            });
        nextButton.append("span")
            .text("Next ");
        nextButton.append("span")
            .attr("class", "glyphicon glyphicon-arrow-right");

      function currentItem() {
          if (selector !== undefined) {
            return list[currentChromIndex][selector];
          }
          else {
            return list[currentChromIndex];
          }
      }

      model.on('change', function() {
          container.select("text").text(currentItem);
      });

      return this;
    }
  });

  return ListSelectorView;

}(d3, utils));
