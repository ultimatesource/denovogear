var store = (function(Immutable) {
  "use strict";
  
  var theStore = Redux.createStore(reducer);

  var StateRecord = Immutable.Record({
    activeNodeSelection: null,
    activeNode: null,
    showSampleTrees: false,
    mutationRecordIndex: 0
  });
  var initialState = new StateRecord();
 
  function reducer(state = initialState, action) {
    switch(action.type) {
      case "NODE_CLICKED":
        return state
          .set('activeNodeSelection', action.selection)
          .set('activeNode', action.node);
      case "TOGGLE_SAMPLE_TREES":
        return state.set('showSampleTrees', !state.showSampleTrees);
      case "MUTATION_CLICKED":
        return state.set('mutationRecordIndex', action.mutationRecordIndex);
      default:
        console.log("Unknown action:", action);
        return state;
    }
  }

  return theStore;
}(Immutable));