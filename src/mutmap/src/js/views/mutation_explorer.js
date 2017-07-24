// eslint exceptions
//
/* global d3 */
/* global PubSub */
/* global pedParser */
/* global pedigreeView*/
/* global vcfParser */
/* global pedigr */
/* global utils */
/* global contigView */

/* global pedigreeFileText */
/* global layoutData */

var mutationExplorerView = (function(d3, PubSub, utils) {
  "use strict";

  function MutationExplorerView(options) {

    utils.optionsManager.checkOptions({
      requiredOptions: ['renderInto', 'graphData', 'pedGraph', 'vcfData'],
      providedOptions: options
    });

    var graphData = options.graphData;
    var pedGraph = options.pedGraph;
    var vcfData = options.vcfData;

    // TODO: Using globals. Hack. Use a better method.
    var selectedContigIndex = 0;
    var selectedMutationIndex = 0;
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

    // Create genome browser view
    var v = contigView.createContigView({
      renderInto: genomeBrowserRenderElement,
      vcfData: vcfData,
      contigData: contigData,
      selectedContigIndex: selectedContigIndex,
      selectedMutationIndex: selectedMutationIndex
    });

    //console.log(vcfData);

    dngOverlay(vcfData.header, vcfData.records[0]);

    var pedigreeRow = parent.append("div")
        .attr("class", "row pedigree-row");

    var pedigreeContainer = pedigreeRow.append("div")
        .attr("class",
              "pedigree-container col-xs-12 col-md-8 panel panel-default"); 

    var pedigreeViewOptions = {
      renderInto: pedigreeContainer,
      graphData: graphData,
    };
    var pedView = pedigreeView.createPedigreeView(pedigreeViewOptions);

    var statsColumn = pedigreeRow.append("div")
        .attr("class", "col-xs-4 col-md-4 stats-col");

    var statsContainer = statsColumn.append("div")
        .attr("class", "stats-container panel panel-default");

    new mutmap.StatsView({
      el: statsContainer
    });

    statsColumn.append("button")
        .attr("id", "sample_tree_toggle")
        .attr("type", "button")
        .attr("class", "btn btn-success")
        .text("Show Trees")
        .on("click", function() {
          PubSub.publish("SAMPLE_TREE_TOGGLE");
        });


    window.addEventListener("resize", function() {
      var dimensions = windowDimensions();

      PubSub.publish("WINDOW_RESIZE");
    });

    PubSub.subscribe("MUTATION_CLICKED", function(topic, data) {
      updateMutation();
    });

    PubSub.subscribe("PREV_MUTATION_BUTTON_CLICKED", function(topic, data) {
      
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

      updateMutation();
    });

    PubSub.subscribe("NEXT_MUTATION_BUTTON_CLICKED", function(topic, data) {

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

      updateMutation();
    });

    function updateMutation() {
      ownerParentageLink.getData().mutation = undefined;
      dngOverlay(vcfData.header,
        contigData[selectedContigIndex].records[selectedMutationIndex]);
      PubSub.publish("MUTATION_INDEX_UPDATED", { 
        selectedContigIndex: selectedContigIndex,
        selectedMutationIndex: selectedMutationIndex
      });
      PubSub.publish("DNG_OVERLAY_UPDATE");
    }

    PubSub.subscribe("MUTATION_SELECTED", function(topic, data) {

      // find matching contig and mutation indexes
      for (var i = 0; i < contigData.length; i++) {
        if (contigData[i].id === data.CHROM) {

          var contigIndex = i;
          for (var j = 0; j < contigData[i].records.length; j++) {
            if (contigData[i].records[j].POS === data.POS) {
              selectedContigIndex = i;
              selectedMutationIndex = j;
              updateMutation();
              break;
            }
          }
          break;
        }
      }
    });

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
 
}(d3, PubSub, utils));
