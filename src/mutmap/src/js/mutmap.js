(function() {

  // The placeholder tags below will be replaced with the correct data objects
  // by the build system. See mutmap/tools/layout_and_template.R
  //
  // provides pedigreeFileText
  /*PEDIGREE_FILE_TEXT_PLACEHOLDER*/
  // provides layoutData
  /*LAYOUT_DATA_PLACEHOLDER*/
  // provides dngOutputFileText
  /*DNG_VCF_DATA_PLACEHOLDER*/

  window.pedigreeFileText = pedigreeFileText;
  window.layoutData = layoutData;
  window.dngOutputFileText = dngOutputFileText;

  //$('.nav-tabs a[href="#mutation-distribution"]')
  //  .on('shown.bs.tab', function(e) {

  var options = {
    renderInto: d3.select(".mutation-distribution-container"),
  };
  mutationDistributionView.createMutationDistributionView(options);

  //});

  var rendered = false;

  $('.nav-tabs a[href="#mutation-explorer"]').on('shown.bs.tab', function(e) {

    if (!rendered) {
      rendered = true;

      var options = {
        renderInto: d3.select(".mutation-explorer-container"),
      };
      mutationExplorerView.createMutationExplorerView(options);
    }

  });


}());
