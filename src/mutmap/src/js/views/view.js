var mutmap = mutmap || {};

mutmap.View = (function(d3, PubSub, utils) {
  "use strict";

  function View(options) {

    utils.optionsManager.checkOptions({
      requiredOptions: ['renderInto'],
      providedOptions: options
    });

    this._renderInto = options.renderInto;
  }

  return View;
  
}(d3, PubSub, utils));
