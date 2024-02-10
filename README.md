# rocktools

Tools for creating and manipulating triangular meshes

Copyright (C) 1999-2017  Mark J. Stock

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


-----------------------------------------------------------------------

## Table of Contents

1. Introduction
2. Getting Started
3. File Summary
4. Rocktools Usage
    1. rockcreate
    2. rockdetail
    3. rocksmooth
    4. rocktrim
    5. rockerode
    6. rockxray
5. File Formats
6. Version History
7. Credits
8. Appendixes


-----------------------------------------------------------------------

## 1.0 Introduction

Rocktools is a package of related programs useful for defining natural
and fractal-like shapes such as those of rocks and rough terrain with
irregular meshes of triangles.

The project was initiated in February of 1999 by Mark Stock to fill
an apparent gap in current, free software to allow recursive roughening
of arbitrary triangle meshes. It is intended to be used for, but not
neccessarily limited to, creation of realistic-looking rock shapes and,
possibly, landscapes.

If you find rocktools useful, or use it in a competition, the author
would appreciate hearing from you.

-----------------------------------------------------------------------

## 2.0 Getting Started

To get rocktools up and running, you only need to build the software
and place the executables in an appropriate place. This may be as
easy as running:

    make

in the directory with this file.

If you downloaded the basic source code distribution, you will need
to build the package. On a Unix system, look through the file called
`Makefile` to see if all of the parameters are set properly. If so,
simply type `make` at the command-line in the src/ directory and the
tools will be built.

To install the files in `/usr/local/bin` (where they will be visible to
a normal user), edit the `Makefile` to expose the "BIN" variable and then

    make
    make install

To learn about the command-line options available and the usage of
any program, just run the program with no command-line options or
with `-help`.


-----------------------------------------------------------------------

## 3.0 File Summary

Rocktools, at the moment, is composed of five programs, summarized below:

* **rockconvert** - Convert between the (few) file formats that Rocktools
    supports.

* **rockcreate** - Create an initial rock shape by wrapping a convex
    hull around randomly created points.

* **rockdetail** - Recursively detail a rock surface composed
    of an irregular triangle mesh.

* **rockdice** - Recursively split a large tri mesh into smaller chunks
    with smooth edges

* **rockerode** - Perform a simple erosion routine over the tri mesh.

* **rockinfo** - dump the node and tri counts, plus min/max of the nodes.

* **rockmarker** - Place simple objects onto a trimesh.

* **rocksmooth** - Smooth a triangle mesh using a Laplace-like operator.

* **rocksplit** - Split a tri mesh into two meshes along a x, y, or z plane

* **rocktrim** - Split an arbitraty tri mesh along a coordinate plane.

* **rockxray** - Write a grayscale image of the shell of a tri mesh.

* **rockpng** - generate a mesh from a heightfield image

* **rockslice** - cut an axis-aligned 2D slice from a trimesh

* **rockbalance** - find the most stable orientation of a closed trimesh

* **rockbob** - create a brick-of-bytes voxel file from a trimesh

Some other C files support the main programs listed above. These files
and descriptions follow:

* `NAMEutil.c` - Includes subroutines used in rockNAME

* `structs.h` - A pseudo-header file describing the dominant data
    structures used in the programs

* `inout.c` - All of the triangle mesh file input and output is handled
    by routines in this file

* `utils.c` - Some vector math and triangle mesh algorithms are kept
    in this file

In addition, there are several files included in the distribution,
these are summarized and described below:

* `README.md` - this documentation file

* `Makefile` - this is a standard Unix Makefile, running "make" on
    a Unix system will use information in this file to build the
    binary executables for Rocktools

Also, there are some basic building-block meshes that serve as good
initial input to rockdetail:

* `onetri.obj`, `square.obj` - one and two triangles, respectively

* `tetra.obj` - a starter triangle mesh file describing a tetrahedron

* `cube.obj` - a starter triangle mesh file describing a cube in 12 triangles

* `rook.obj` - a Wavefront file of the chesspiece of the same name

* `icosahedron.obj` - a 20-sided solid, useful for recursively
    refining a sphere using rockdetail's `-sph` option


-----------------------------------------------------------------------

## 4.0 Program Usage

All programs accept these two arguments:

    -okey       specify output format, key= raw, rad, pov, obj, tin, rib    
                default = raw; surface normal vectors are not supported     
                                                                           
    -h
    -help       (anywhere in the command) returns this help information     
 
Also, all options may be abbreviated to their unambiguous length,
and all triangle mesh output is to stdout. stderr contains status
messages that, as of this writing, cannot be turned off.

Programs that take an input file can accept triangle mesh files with the
following formats: .obj, .raw, or .tin. The programs require the input
file to use its valid 3-character filename extension.

To see the usage for each program, run the executable with no options, or
with the `-help` option.


-----------------------------------------------------------------------

## 5.0 File Formats

Rocktools can read three and write six different file formats. They are
described below. Keep in mind that in the actual files, the `x1 y1 z1`
notations would be replaced with actual floating-point numbers.

* Raw Triangle Format (.raw) - the most basic triangle mesh file
    format possible, it consists of nine numbers, corresponding
    to three nodes of a triangle in three dimensions, per line.
    Comments are supported with a '#' as the first non-whitespace
    character on a line.

    The .tin and .raw formats are currently the only supported formats
    for input into rockdetail/rocktools.

        x1 y1 z1 x2 y2 z2 x3 y3 z3

* Triangle Irregular Mesh (.tin) - a file format used by scape-1.2, a
    surface simplification tool written by Michael Garland, found at:
    http://www.cs.cmu.edu/afs/cs/user/garland/www/scape/index.html

    The .tin and .raw formats are currently the only supported formats
    for input into rockdetail/rocktools.

    Each triangle is entered on its own line, in the format

        t x1 y1 z1 x2 y2 z2 x3 y3 z3

* Wavefront (.obj) - a standard 3D file format used by Alias Wavefront,
    it is a free-form file which lists coordinates of all vertexes and
    each face references the vertex numbers. An example follows:

        v x1 y1 z1
        v x2 y2 z2
        v x3 y3 z3
        f 1 2 3

    This is the recommended file format, as the reading routines are
    optimized for very large files (1/2 million tris+).

* Persistence of Vision Raytracer (.pov) - a file format used by the 
    popular raytracer, POV, the file format is very similar to .tin:

        mesh {
        triangle {<x1 y1 z1> <x2 y2 z2> <x3 y3 z3>}
        }

* Radiance (.rad) - a scene file format used by Radiance, a very accurate
    lighting simulation software package. Also similar to .tin. Note that
    the inherited surface attribute is called "default."
    http://radsite.lbl.gov/radiance/HOME.html

        default polygon p0
        0
        0
        9 x1 y1 z1 x2 y2 z2 x3 y3 z3

* Renderman (.rib) - A Renderman-compliant polygon definition, used by
    the Blue Moon RayTracer (BMRT) and Pixar/Gritz's Renderman. It
    has one triangle per line and goes a little like this:

        Polygon "P" [ x1 y1 z1 x2 y2 z2 x3 y3 z3 ]


-----------------------------------------------------------------------

## 6.0 Version History

    1999-02-24  v0.0    Started project
    1999-03-11  v0.1    rockdetail works, outputs POV, Rad, Wavefront, .tin
    1999-03-20  v0.1.1  added RIB and RAW output, reworked input/output
    1999-03-31  v0.1.2  added clamp_edges option to rockdetail
    1999-04-20  v0.2a   rockcreate added
    1999-05-03  v0.3a   rocksmooth added, reorganized directory
    1999-05-28  v0.4a   rocktrim added, new recursion options in rockdetail
    2000-12-02  v0.5a   rockerode added
    2003-09-16  v0.6    obj input, new options to rockdetail
    2004-01-26  v0.7    faster obj input, added rockxray, minor other fixes
    2006-06-26  v0.8    added more features, added rockmarker
    2007-11-23  v0.9    small fixes, reduced memory use, broke rockcreate
    2009-04-22  v1.0    rewrote rockcreate with new features, added rockinfo,
                        rocksplit, and rockdice
    2011-08-01  v1.1    added rockpng, to convert 16-bit PNGs to meshes
    2014-04-06  v1.2    added rockslice, new features to rockpng
    2017-09-12  v1.3    many improvements to xray, moved to github

-----------------------------------------------------------------------

## 7.0 Credits

My thanks to all those who helped or showed interest in rocktools.
Their enthusiasm kept me interested too, and cannot be thanked enough.
The GNU Public License protects their work in rocktools, too.

* Jeff Balcerski - Testing, POV support, images/movies
* Anthony D'Agostino - DOS port, testing, code fixes
* Steven Pigeon - Initial concept, page design
* Ken Clarkson - "hull" program
* Jeff Norris - "Alpha" testing

## Citing Rocktools

I don't get paid for writing or maintaining this, so if you find this tool useful or mention it in your writing, please please cite it by using the following BibTeX entry.

```
@Misc{Rocktools2017,
  author =       {Mark J.~Stock},
  title =        {Rocktools:  Tools for creating and manipulating triangle meshes},
  howpublished = {\url{https://github.com/markstock/rocktools}},
  year =         {2017}
}
```
