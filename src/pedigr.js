(function() {
  "use strict";

  function PedigreeGraph() {
    this.nodes = {};
  };

  PedigreeGraph.prototype.getPerson = function(id) {
    return this.nodes[id];
  };

  PedigreeGraph.prototype.getPersons = function(query) {
    var persons = [];
    for (var key in this.nodes) {
      persons.push(this.nodes[key]);
    }
    return persons;
  };

  PedigreeGraph.prototype.addPerson = function(attr) {
    this.nodes[attr.id] = Person.createPerson(attr);
    return this.nodes[attr.id];
  };

  PedigreeGraph.createGraph = function() {
    var graph = new PedigreeGraph();
    return graph;
  };

  function Person(attr) {
    this.id = attr.id;
    this.sex = attr.sex;
    this.father = attr.father;
    this.mother = attr.mother;
    this.marriage = attr.marriage;
    this.data = attr.data;
  };

  Person.createPerson = function(attr) {
    return new Person(attr);
  };

  function PersonBuilder(id) {
    this._id = id;
    this._sex = undefined;
    this._father = undefined;
    this._mother = undefined;
    this._marriage = undefined;
    this._data = undefined;
    return this;
  };

  PersonBuilder.prototype.sex = function(sex) {
    this._sex = sex;
    return this;
  }

  PersonBuilder.prototype.father = function(father) {
    this._father = father;
    return this;
  };

  PersonBuilder.prototype.mother = function(mother) {
    this._mother = mother;
    return this;
  };

  PersonBuilder.prototype.marriage = function(marriage) {
    this._marriage = marriage;
    return this;
  };

  PersonBuilder.prototype.data = function(data) {
    this._data = data;
    return this;
  };

  PersonBuilder.prototype.build = function() {
    return Person.createPerson({
      id: this._id,
      sex: this._sex,
      father: this._father,
      mother: this._mother,
      marriage: this._marriage,
      data: this._data
    });
  };

  PersonBuilder.createPersonBuilder = function(id) {
    return new PersonBuilder(id);
  };


  module.exports = {
    PedigreeGraph: PedigreeGraph,
    PersonBuilder: PersonBuilder
  };
}());
