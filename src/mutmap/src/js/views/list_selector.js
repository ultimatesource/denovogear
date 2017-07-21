mutmap.ListSelectorView = (function(d3, utils) {
  "use strict";
  
  var ListSelectorView = Backbone.View.extend({
    
    initialize: function(options) {

      this.list = options.list;
      this.selector = options.selector;
      this.itemName = options.itemName;

      this.currentChromIndex = 0;

      this.container = options.el.append("div")
          .attr("class", "list-selector");

      this.model.on('change', function() {
          this.container.select("text").text(this.currentItem.bind(this));
      }, this);

      this.render();
    },

    render: function() {

      // Preserve this context for use in anonymous functions
      var self = this;

      this.container.append("div")
        .append("text")
          .text(this.currentItem.bind(this));

      var prevButton = this.container.append("button")
          .attr("class", "btn btn-default")
          .on("click", function() {
            self.currentChromIndex--;
            if (self.currentChromIndex < 0) {
              self.currentChromIndex = self.list.length - 1;
            }

            self.model.set('index', self.currentChromIndex);
          });
      prevButton.append("span")
          .attr("class", "glyphicon glyphicon-arrow-left");
      prevButton.append("span")
          .text(" Previous");

      var dropdown = this.container.append("div")
          .attr("class", "dropdown")
      dropdown.append("button")
          .attr("class", "btn btn-default dropdown-toggle")
          .attr("id", "dropdownMenu1")
          .attr("type", "button")
          .attr("data-toggle", "dropdown")
          .text("Select " + this.itemName);

      dropdown.append("ul")
          .attr("class", "dropdown-menu")
        .selectAll(".mutation")
          .data(this.list)
        .enter().append("li")
          .attr("class", "mutation")
        .append("a")
          .attr("href", "#")
          .text(function(d) {
            if (self.selector !== undefined) {
              return d[self.selector];
            }
            else {
              return d;
            }
          })
          .on("click", function(d, i) {
            self.currentChromIndex = i;
            self.model.set('index', self.currentChromIndex);
          });

        var nextButton = this.container.append("button")
            .attr("class", "btn btn-default")
            .on("click", function() {
              self.currentChromIndex++;
              if (self.currentChromIndex === self.list.length) {
                self.currentChromIndex = 0;
              }

              self.model.set('index', self.currentChromIndex);
            });
        nextButton.append("span")
            .text("Next ");
        nextButton.append("span")
            .attr("class", "glyphicon glyphicon-arrow-right");


      return this;
    },

    currentItem: function() {
      if (this.selector !== undefined) {
        return this.list[this.currentChromIndex][this.selector];
      }
      else {
        return this.list[this.currentChromIndex];
      }
    }
  });

  return ListSelectorView;

}(d3, utils));
