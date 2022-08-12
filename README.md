# icosphere
Creating a Class I geodesic icosahedron of unit radius through subdivision.

This MATLAB script contains options to choose two different starting regular icosahedron at different rotation angles, and two different subdivision algorithms.

## Usage

`[V,F] = icosphere(N);` generates to matrices containing vertices `V` and faces `F` in the form of rowID in `V`, with each edge subdivided into `N` equal segments.

Use `trimesh(icosphere(N));axis equal;` to plot the resulting sphere.

## Geodesic Notation
This script generates a $$\{3,5+\}_{(n,0)}$$ [geodesic polyhedron](https://en.wikipedia.org/wiki/Geodesic_polyhedron).

## References
This script is inspired by 

- [C# code by Andres Kahler](http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html)
- [MATLAB code by Wil O.C. Ward](https://uk.mathworks.com/matlabcentral/fileexchange/50105-icosphere)
- [Python Geodesic Icosphere package by vedranaa](https://github.com/vedranaa/icosphere)
