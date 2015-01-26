![Unfolded sphere](https://farm8.staticflickr.com/7464/15880804609_88b0865685_z.jpg)

`unfold` reads a binary STL file on standard input and generates a
SVG that contains the triangles "folded flat" so that they can be
laser cut.  It will output multiple groups in the SVG file that
will need to be re-arranged to fit on the laser cutter bed.

More info: https://trmm.net/Unfolding_STL

This is a work in progress -- it is not yet feature complete. Current features:

* Collision detection works to build a non-overlapping set of triangles.

* Coplanar triangles will be merged if possible.

* Mountain/valley folds are marked.  Dashed line means valley fold.

* Starting face can be selected or randomly chosen; some produce better results than others.

* `stl-convert` script can convert OpenSCAD ASCII STL files into binary STL files for `unfold` to process.

Among the features that it could use:

* A better heuristic for finding the maximum non-overlaping set of triangles
(Currently breadth-first search is used, with a slight preference for coplanar
triangles)

* Tabs for securing parts together.

* Collapsing of very small or very thin triangles.

![Wireframed dodecahedron](https://farm8.staticflickr.com/7320/16181756310_729bb6e186_z.jpg)

`wireframe` reads a binary STL file on standard input and generates an
OpenSCAD file that has each vertex replaced with a 3D printable
connector to create a large-scale version of the object.

More info: https://trmm.net/Wireframe

This is very beta! It desperately needs some command line options and
some smarts in labeling the connectors, especially if the structure
is not regular.  Also, it does not detect intersections between
beams or if the nodes are too close together.
