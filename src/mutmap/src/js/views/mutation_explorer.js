var mutmap = mutmap || {};

(function() {
  "use strict";

  mutmap.MutationExplorerView = Backbone.View.extend({
    initialize: function(options) {

      this.pedGraph = options.pedGraph;
      this.vcfData = options.vcfData;

      this.ownerParentageLink;

      // transform contigs to hierarchical format
      this.contigData = [];
      this.vcfData.header.contig.forEach(function(contig) {
        var id = contig.ID;
        var contig = {
          id: id,
          length: contig.length,
          records: []
        };

        this.vcfData.records.forEach(function(record) {
          if (id === record.CHROM) {
            contig.records.push(record);
          }
        });

        this.contigData.push(contig);
      }, this);

      this.d3el = d3.select(this.el);

      var dim = utils.getDimensions(this.d3el);

      var genomeBrowserRenderElement = this.d3el.append("div")
          .attr("class", "row")
        .append("div")
          .attr("id", "genome_browser_wrapper")
          .attr("class", "genome-browser-container col-xs-12 panel panel-default");

      var contigDimensions = utils.getDimensions(genomeBrowserRenderElement);

      var ContigModel = Backbone.Model.extend({
        defaults: {
          selectedContigIndex: 0,
          selectedMutationIndex: 0,
          contigData: this.contigData[0]
        }
      });

      this.contigModel = new ContigModel();

      // Create genome browser view
      new mutmap.ContigView({
        el: genomeBrowserRenderElement.node(),
        model: this.contigModel,
        prevMutationButtonClicked: this.prevMutationButtonClicked.bind(this),
        nextMutationButtonClicked: this.nextMutationButtonClicked.bind(this),
        mutationClicked: this.mutationClicked.bind(this)
      });

      this.dngOverlay(this.vcfData.records[0]);

      var ActiveNodeModel = Backbone.Model.extend({
        defaults: {
          activeNode: null,
        }
      });

      this.activeNodeModel = new ActiveNodeModel();

      var pedigreeRow = this.d3el.append("div")
          .attr("class", "row pedigree-row");

      var pedigreeContainer = pedigreeRow.append("div")
          .attr("class", "col-xs-12 col-md-8 panel panel-default")
            .style("height", (dim.height - contigDimensions.height)+'px');

      var PedigreeModel = Backbone.Model.extend({
        defaults: {
          activeNodeModel: this.activeNodeModel,
          graphData: options.graphData,
          showSampleTrees: false
        }
      });

      this.pedigreeModel = new PedigreeModel();

      new mutmap.PedigreeView({
        el: pedigreeContainer,
        model: this.pedigreeModel,
      });

      var statsColumn = pedigreeRow.append("div")
          .attr("class", "col-xs-4 col-md-4 stats-col");

      var statsContainer = statsColumn.append("div")
          .attr("class", "stats-container panel panel-default");

      new mutmap.StatsView({
        el: statsContainer,
        model: this.activeNodeModel,
      });

      var context = this;
      statsColumn.append("button")
          .attr("id", "sample_tree_toggle")
          .attr("type", "button")
          .attr("class", "btn btn-success")
          .text("Show Trees")
          .on("click", function() {
            context.pedigreeModel.set('showSampleTrees',
              !context.pedigreeModel.get('showSampleTrees'));
          });
    },

    prevMutationButtonClicked: function() {

      var selectedContigIndex = this.contigModel.get('selectedContigIndex');
      var selectedMutationIndex = this.contigModel.get('selectedMutationIndex');
      
      selectedMutationIndex--;
      if (selectedMutationIndex === -1) {

        selectedContigIndex--;
        if (selectedContigIndex === -1) {
          selectedContigIndex = this.contigData.length - 1;
          selectedMutationIndex = 0;
        }

        selectedMutationIndex =
          this.contigData[selectedContigIndex].records.length - 1;
      }

      this.contigModel.set('selectedContigIndex', selectedContigIndex);
      this.contigModel.set('selectedMutationIndex', selectedMutationIndex);
      this.contigModel.set('contigData', this.contigData[selectedContigIndex]);

      this.updateMutation();
    },

    nextMutationButtonClicked: function() {

      var selectedContigIndex = this.contigModel.get('selectedContigIndex');
      var selectedMutationIndex =
        this.contigModel.get('selectedMutationIndex');

      selectedMutationIndex++;
      if (selectedMutationIndex ===
          this.contigData[selectedContigIndex].records.length) {

        selectedMutationIndex = 0;

        selectedContigIndex++;
        if (selectedContigIndex === this.contigData.length) {
          selectedContigIndex = 0;
          selectedMutationIndex = 0;
        }
      }

      this.contigModel.set('selectedContigIndex', selectedContigIndex);
      this.contigModel.set('selectedMutationIndex', selectedMutationIndex);
      this.contigModel.set('contigData', this.contigData[selectedContigIndex]);

      this.updateMutation();
    },

    updateMutation: function() {
      var selectedContigIndex = this.contigModel.get('selectedContigIndex');
      var selectedMutationIndex =
        this.contigModel.get('selectedMutationIndex');

      this.ownerParentageLink.getData().mutation = undefined;
      this.dngOverlay(
        this.contigData[selectedContigIndex].records[selectedMutationIndex]);

      // TODO: This feels like a hack. At least having both of them does.
      this.pedigreeModel.trigger('change');
      this.activeNodeModel.trigger('change');
    },

    mutationClicked: function(data) {

      var selectedContigIndex = this.contigModel.get('selectedContigIndex');
      var selectedMutationIndex =
        this.contigModel.get('selectedMutationIndex');

      // find matching contig and mutation indexes
      for (var i = 0; i < this.contigData.length; i++) {
        if (this.contigData[i].id === data.CHROM) {

          var contigIndex = i;
          for (var j = 0; j < this.contigData[i].records.length; j++) {
            if (this.contigData[i].records[j].POS === data.POS) {
              selectedContigIndex = i;
              selectedMutationIndex = j;

              this.contigModel.set('selectedContigIndex', selectedContigIndex);
              this.contigModel.set('selectedMutationIndex',
                selectedMutationIndex);
              this.updateMutation();
              break;
            }
          }
          break;
        }
      }
    },

    dngOverlay: function(record) {

      var mutationLocation = record.INFO.DNL.slice(3);
      var owner = findOwnerNode(this.pedGraph, mutationLocation);

      if (owner !== undefined) {
        this.ownerParentageLink = owner.getParentageLink();
        var parentageLinkData = {
          mutation: record.INFO.DNT
        };
        this.ownerParentageLink.setData(parentageLinkData);

        this.vcfData.header.sampleNames.forEach(function(sampleName) {
          var format = record[sampleName];


          // TODO: Using libraries for now, but might be more correct to use
          // GL-1, GL-2, etc nodes?
          if (isLibraryNode(sampleName)) {

            var id = getIdFromLibraryName(sampleName);
            var personNode = this.pedGraph.getPerson(id);
            personNode.data.dngOutputData = format;
          }

        }, this);
      }
      else {
        throw "No mutation found";
      }
    },

  });

  function findOwnerNode(pedGraph, sampleName) {
    //var strippedName = getStrippedName(sampleName);
    var persons = pedGraph.getPersons();
    for (var index = 0; index < persons.length; index++) {
      var person = persons[index];
      //var sampleNode = findInTree(person.data.sampleIds, strippedName);
      var sampleNode = findInTree(person.data.sampleIds, sampleName);
      if (sampleNode !== undefined) {
        return person;
      }
    }
    return undefined;
  }

  function findInTree(tree, sampleName) {

    // TODO: This seems very likely to break in the future. Need to find a
    // robust way of matching up sample names and libraries.
    if (tree.name != "" && tree.name.includes(sampleName)) {
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

  function isLibraryNode(sampleName) {
    return sampleName.startsWith("LB/");
  }

  function getIdFromSampleName(sampleName) {
    return sampleName.slice(3);
  }

  function getIdFromLibraryName(sampleName) {
    //return Number(sampleName.slice(-3));
    return sampleName.slice(3, 10);
  }

  function createMutationExplorerView(options) {
    return new MutationExplorerView(options);
  }

  return {
    createMutationExplorerView: createMutationExplorerView
  };
 
}());
