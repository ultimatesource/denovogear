(function() {
  "use strict";

  function PedigreeGraph() {
    this.persons = {};
    this.marriages = [];
    this.parentageLinks = {};
  };

  PedigreeGraph.prototype.getPerson = function(id) {
    return this.persons[id];
  };

  PedigreeGraph.prototype.getPersons = function(query) {
    var persons = [];
    for (var key in this.persons) {
      persons.push(this.persons[key]);
    }
    return persons;
  };

  PedigreeGraph.prototype.addPerson = function(person) {
    this.persons[person.id] = person;
    return person;
  };

  PedigreeGraph.prototype.addMarriage = function(marriage) {
    this.marriages.push(marriage);
    return marriage;
  };

  PedigreeGraph.createGraph = function() {
    var graph = new PedigreeGraph();
    return graph;
  };

  function Person(attr) {
    this.id = attr.id;
    this.sex = attr.sex;
    //this.father = attr.father;
    //this.mother = attr.mother;
    this.parentageLink = attr.parentageLink;
    this.marriageLink = attr.marriageLink;
    this.data = attr.data;
  };

  Person.createPerson = function(attr) {
    return new Person(attr);
  };

  Person.prototype.getParentageLink = function() {
    return this.parentageLink;
  };

  function PersonBuilder(id) {
    this._id = id;
    this._sex = undefined;
    //this._father = undefined;
    //this._mother = undefined;
    this._parentageLink = undefined;
    this._marriageLink = undefined;
    this._data = undefined;
    return this;
  };

  PersonBuilder.prototype.sex = function(sex) {
    this._sex = sex;
    return this;
  }

  PersonBuilder.prototype.parentageLink = function(parentageLink) {
    this._parentageLink = parentageLink;
    return this;
  };

  //PersonBuilder.prototype.father = function(father) {
  //  this._father = father;
  //  return this;
  //};

  //PersonBuilder.prototype.mother = function(mother) {
  //  this._mother = mother;
  //  return this;
  //};

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
      //father: this._father,
      //mother: this._mother,
      parentageLink: this._parentageLink,
      marriageLink: this._marriageLink,
      data: this._data
    });
  };

  PersonBuilder.createPersonBuilder = function(id) {
    return new PersonBuilder(id);
  };

  function Marriage(attr) {
    this.fatherLink = attr.fatherLink;
    this.motherLink = attr.motherLink;
    this.childLinks = attr.childLinks;
    this.data = attr.data;
  }

  Marriage.prototype.addSpouse = function(spouse) {
    if (spouse.sex === 'male') {
      var link = MarriageLink.createMarriageLink(this, spouse);
      this.fatherLink = link;
      spouse.marriageLink = link;
    }
    else {
      var link = MarriageLink.createMarriageLink(this, spouse);
      this.motherLink =  link;
      spouse.marriageLink = link;
    }
  };

  Marriage.prototype.addChild = function(child) {
    var parentageLink = ParentageLink.createParentageLink(this, child);
    this.childLinks.push(parentageLink);
    child.parentageLink = parentageLink;
    return parentageLink;
  };

  Marriage.createMarriage = function(attr) {
    return new Marriage(attr);
  };

  function MarriageBuilder() {
    this._father = undefined;
    this._mother = undefined;
    this._children = undefined;
    this._data = undefined;
  }

  MarriageBuilder.createMarriageBuilder = function() {
    return new MarriageBuilder();
  };

  MarriageBuilder.prototype.spouse = function(spouse) {
    if (spouse.sex === 'male') {
      this._father = spouse;
    }
    else {
      this._mother = spouse;
    }
    return this;
  };

  MarriageBuilder.prototype.father = function(father) {
    this._fatherLink = father;
    return this;
  };

  MarriageBuilder.prototype.mother = function(mother) {
    this._mother = mother;
    return this;
  };

  MarriageBuilder.prototype.children = function(children) {
    this._children = children;
    return this;
  };

  MarriageBuilder.prototype.data = function(data) {
    this._data = data;
    return this;
  };

  MarriageBuilder.prototype.build = function() {
    var marriage = new Marriage.createMarriage({
      data: this._data
    });

    marriage.addSpouse(this._father);
    marriage.addSpouse(this._mother);
    marriage.childLinks = [];
    return marriage;
  }

  function MarriageLink(marriage, spouse) {
    this.marriage = marriage;
    this.spouse = spouse;
  }

  MarriageLink.createMarriageLink = function(marriage, spouse) {
    return new MarriageLink(marriage, spouse);
  };

  function ParentageLink(parentage, child) {
    this.parentage = parentage;
    this.child = child;
    this.data = undefined;
  }

  ParentageLink.createParentageLink = function(parentage, child) {
    return new ParentageLink(parentage, child);
  };

  ParentageLink.prototype.setData = function(data) {
    this.data = data;
  };


  module.exports = {
    PedigreeGraph: PedigreeGraph,
    PersonBuilder: PersonBuilder,
    MarriageBuilder: MarriageBuilder
  };

}());
