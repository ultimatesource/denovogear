(function() {
  "use strict";

  function PedigreeGraph() {
    this.nodes = {};
  }

  PedigreeGraph.prototype.getPerson = function(id) {
    return this.nodes[id];
  }

  PedigreeGraph.prototype.addPerson = function(attr) {
    this.nodes[attr.id] = Person.createPerson(attr);
    return this.nodes[attr.id];
  }

  PedigreeGraph.createGraph = function() {
    var graph = new PedigreeGraph();
    return graph;
  }

  function Person(attr) {
    this.id = attr.id;
    this.sex = attr.sex;
    this.father = attr.father;
    this.mother = attr.mother;
    this.marriage = attr.marriage;
    this.data = attr.data;
  }

  Person.createPerson = function(attr) {
    return new Person(attr);
  }


  module.exports = {
    PedigreeGraph: PedigreeGraph,
    Person: Person
  };
}());
