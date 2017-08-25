# Docs

# View width and height

View dimensions are computed by the owner of the view object (often a parent
view). Once computed a container of the proper dimensions is created by the
owner, and passed in as `el`. The convention is for each view to fill the
provided element.

HTML elements (usually `div` elements) use CSS styling to control the width and
height.  SVG elements set the width and height directly. Since SVG `g` elements
cannot directly have width and height attributes, nested `svg` elements are
used.
