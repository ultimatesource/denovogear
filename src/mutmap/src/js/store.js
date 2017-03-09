var store = (function() {
  
  var theStore = Redux.createStore(reducer);

  var StateRecord = Immutable.Record({
    activeNodeSelection: null,
    activeNode: null,
    showSampleTrees: false
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
        console.log(action);
        //return state.set();
        return state;
      default:
        console.log("Invalid action");
        return state;
    }
  }

  return theStore;
}());
