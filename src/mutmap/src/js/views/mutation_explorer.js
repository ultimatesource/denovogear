// eslint exceptions
//
/* global d3 */
/* global pedParser */
/* global pedigreeView*/
/* global vcfParser */
/* global pedigr */
/* global utils */
/* global contigView */

/* global pedigreeFileText */
/* global layoutData */

var mutationExplorerView = (function(d3, utils) {
  "use strict";

  function MutationExplorerView(options) {

    utils.optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'graphData', 'pedGraph', 'vcfData'],
      providedOptions: options
    });

    var pedGraph = options.pedGraph;
    var vcfData = options.vcfData;

    // TODO: Using globals. Hack. Use a better method.
    var ownerParentageLink;

    // transform contigs to hierarchical format
    var contigData = [];
    vcfData.header.contig.forEach(function(contig) {
      var id = contig.ID;
      var contig = {
        id: id,
        length: contig.length,
        records: []
      };

      vcfData.records.forEach(function(record) {
        if (id === record.CHROM) {
          contig.records.push(record);
        }
      });

      contigData.push(contig);
    });

    var parent = options.renderInto;
    var genomeBrowserRenderElement = parent.append("div")
        .attr("class", "row")
      .append("div")
        .attr("id", "genome_browser_wrapper")
        .attr("class", "genome-browser-container");

    var ContigModel = Backbone.Model.extend({
      defaults: {
        selectedContigIndex: 0,
        selectedMutationIndex: 0,
        contigData: contigData[0]
      }
    });

    var contigModel = new ContigModel();

    // Create genome browser view
    new mutmap.ContigView({
      el: genomeBrowserRenderElement.node(),
      model: contigModel,
      prevMutationButtonClicked: prevMutationButtonClicked,
      nextMutationButtonClicked: nextMutationButtonClicked,
      mutationClicked: mutationClicked
    });

    //console.log(vcfData);

    dngOverlay(vcfData.header, vcfData.records[0]);

    var ActiveNodeModel = Backbone.Model.extend({
      defaults: {
        activeNode: null,
      }
    });

    var activeNodeModel = new ActiveNodeModel();

    var pedigreeRow = parent.append("div")
        .attr("class", "row pedigree-row");

    var pedigreeContainer = pedigreeRow.append("div")
        .attr("class",
              "pedigree-container col-xs-12 col-md-8 panel panel-default"); 

    var PedigreeModel = Backbone.Model.extend({
      defaults: {
        activeNodeModel: activeNodeModel,
        graphData: options.graphData,
        showSampleTrees: false
      }
    });

    var pedigreeModel = new PedigreeModel();

    new mutmap.PedigreeView({
      el: pedigreeContainer,
      model: pedigreeModel,
    });

    var statsColumn = pedigreeRow.append("div")
        .attr("class", "col-xs-4 col-md-4 stats-col");

    var statsContainer = statsColumn.append("div")
        .attr("class", "stats-container panel panel-default");

    new mutmap.StatsView({
      el: statsContainer,
      model: activeNodeModel,
    });

    statsColumn.append("button")
        .attr("id", "sample_tree_toggle")
        .attr("type", "button")
        .attr("class", "btn btn-success")
        .text("Show Trees")
        .on("click", function() {
          pedigreeModel.set('showSampleTrees',
            !pedigreeModel.get('showSampleTrees'));
        });

    function prevMutationButtonClicked() {

      var selectedContigIndex = contigModel.get('selectedContigIndex');
      var selectedMutationIndex = contigModel.get('selectedMutationIndex');
      
      selectedMutationIndex--;
      if (selectedMutationIndex === -1) {

        selectedContigIndex--;
        if (selectedContigIndex === -1) {
          selectedContigIndex = contigData.length - 1;
          selectedMutationIndex = 0;
        }

        selectedMutationIndex =
          contigData[selectedContigIndex].records.length - 1;
      }

      contigModel.set('selectedContigIndex', selectedContigIndex);
      contigModel.set('selectedMutationIndex', selectedMutationIndex);
      contigModel.set('contigData', contigData[selectedContigIndex]);

      updateMutation();
    }

    function nextMutationButtonClicked() {

      var selectedContigIndex = contigModel.get('selectedContigIndex');
      var selectedMutationIndex = contigModel.get('selectedMutationIndex');

      selectedMutationIndex++;
      if (selectedMutationIndex ===
          contigData[selectedContigIndex].records.length) {

        selectedMutationIndex = 0;

        selectedContigIndex++;
        if (selectedContigIndex === contigData.length) {
          selectedContigIndex = 0;
          selectedMutationIndex = 0;
        }
      }

      contigModel.set('selectedContigIndex', selectedContigIndex);
      contigModel.set('selectedMutationIndex', selectedMutationIndex);
      contigModel.set('contigData', contigData[selectedContigIndex]);

      updateMutation();
    }

    function updateMutation() {
      var selectedContigIndex = contigModel.get('selectedContigIndex');
      var selectedMutationIndex = contigModel.get('selectedMutationIndex');

      ownerParentageLink.getData().mutation = undefined;
      dngOverlay(vcfData.header,
        contigData[selectedContigIndex].records[selectedMutationIndex]);

      // TODO: This feels like a hack. At least having both of them does.
      pedigreeModel.trigger('change');
      activeNodeModel.trigger('change');
    }

    function mutationClicked(data) {

      var selectedContigIndex = contigModel.get('selectedContigIndex');
      var selectedMutationIndex = contigModel.get('selectedMutationIndex');

      // find matching contig and mutation indexes
      for (var i = 0; i < contigData.length; i++) {
        if (contigData[i].id === data.CHROM) {

          var contigIndex = i;
          for (var j = 0; j < contigData[i].records.length; j++) {
            if (contigData[i].records[j].POS === data.POS) {
              selectedContigIndex = i;
              selectedMutationIndex = j;

              contigModel.set('selectedContigIndex', selectedContigIndex);
              contigModel.set('selectedMutationIndex', selectedMutationIndex);
              updateMutation();
              break;
            }
          }
          break;
        }
      }
    }

    function dngOverlay(header, record) {

      var mutationLocation = record.INFO.DNL.slice(3);
      var owner = findOwnerNode(mutationLocation);

      if (owner !== undefined) {
        ownerParentageLink = owner.getParentageLink();
        var parentageLinkData = {
          mutation: record.INFO.DNT
        };
        ownerParentageLink.setData(parentageLinkData);

        header.sampleNames.forEach(function(sampleName) {
          var format = record[sampleName];


          // TODO: Using libraries for now, but might be more correct to use
          // GL-1, GL-2, etc nodes?
          if (isLibraryNode(sampleName)) {

            var id = getIdFromLibraryName(sampleName);
            var personNode = pedGraph.getPerson(id);
            personNode.data.dngOutputData = format;
          }

          //if (isPersonNode(sampleName)) {
          //  var id = getIdFromSampleName(sampleName);
          //  var personNode = pedGraph.getPerson(id);
          //  personNode.data.dngOutputData = format;
          //}
          //else {
          //  var sampleNode = findMatchingSampleNode(sampleName);
          //  sampleNode.dngOutputData = format;
          //}
        });
      }
      else {
        throw "No mutation found";
      }
    }

    // TODO this and findMatchingSampleNode have almost the same logic. Find a
    // way to extract the duplication
    function findOwnerNode(sampleName) {
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

    function windowDimensions() {
      var w = window,
          d = document,
          e = d.documentElement,
          g = d.getElementsByTagName('body')[0],
          x = w.innerWidth || e.clientWidth || g.clientWidth,
          y = w.innerHeight|| e.clientHeight|| g.clientHeight;

      return { x: x, y: y };
    }

  }

  function createMutationExplorerView(options) {
    return new MutationExplorerView(options);
  }

  return {
    createMutationExplorerView: createMutationExplorerView
  };
 
}(d3, utils));
