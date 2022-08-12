function [vv,ff] = icosphere(varargin)
%ICOSPHERE Generate a unit geodesic icosahedron.
% Creating a Class I geodesic icosahedron of unit radius through subdivision.
% In geodesic notation, this script generates a {3,5+}_(n,0) geodesic polyhedron.
% (see https://en.wikipedia.org/wiki/Geodesic_polyhedron for more details.)
%
% This script uses two different methods for subvision.
%    If N is a power of 2 (i.e. N=2^m, where m is +ve integer), 
%       the recursive subdivision algorithm is used.
%    Otherwise,
%       the frequency subdivision algorithm is used.
% This script uses two different regular icosahedrons.
%    By default (rotation_flag = false),
%       the regular icosaderhon have vertices x,y,z=0 OR +-1 OR +-(1+sqrt(5)/2).
%    When flagged (rotation_flag = true),
%       the regular icosahedron is rotated with two of the vertices at [+-1,0,0].
%
% Usage:
%   [V,F] = icosphere(N)
%   generates to matrices containing vertices (V)
%   and faces (F) in the form of rowID in V,
%   with each edge subdivided into N equal segments.
%
% Options:
%   
%   TR = ICOSPHERE(N) generates a MATLAB triangulation object.
%
%   [V,F] = icosphere(N, rotation_flag = true)
%   generates the geodesic icosahedron with two of the vertices at [-1,0,0]
%   and [1,0,0].
%
%   [V,F] = icosphere(N, rotation_flag, frequency_flag)
%   enforce algorithm used to subdivision.
%
%   Parameters:
%       N - Number of division on an edge of the regular icosahedron.
%       rotation_flag - false: regular icosahedren; (default)
%                        true: rotated icosahedren with vertex at [1,0,0]
%       freq_div_flag - false: use recursive subdivision algorithm
%                        true: use frequency subdivision algorithm
%                     default: use recursive only when N=2^m
%
%   Inspiration from:
%   C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
%
%   MATLAB code by Wil O.C. Ward 19/03/2015 (U of Nottingham)
%   https://uk.mathworks.com/matlabcentral/fileexchange/50105-icosphere
%
%   Python Geodesic Icosphere package by vedranaa (vand at dtu.dk, 8.12.2017.)
%   https://github.com/vedranaa/icosphere
%   https://doi.org/10.1109/3DV.2014.60


% Parse possible axes input
if nargin > 3
    error('Too many input variables.');
end

n = 3; % default number of sub-divisions
if nargin > 0
    n = double(varargin{1}); % override based on input
    % Error check
    if mod(n,1) || n<=0
        error('N is not an integer');
    end
end 

% generate regular unit icosahedron (20 faced polyhedron)
if nargin >1
    if varargin{2}
        % generate regular unit icosahedron
        % rotated with vertices at [1,0,0] and [-1,0,0]
        [v,f] = icosahedron();
    else
        % generate regulat unit icosahedron,
        [v,f] = icosahedron_alt();
    end
else
    % default
    [v,f] = icosahedron_alt();
end

if nargin > 2
    if varargin{3}
        [v,f]=intdivision(n,f,v);
    else
        [v,f]=subdivision(log2(n),f,v);
    end
else
    % default
    if mod(log2(n),1)
        % If n =/= 2^m
        [v,f]=intdivision(n,f,v);
    else
        % If n = 2^m
        [v,f]=subdivision(log2(n),f,v);
    end
end

switch(nargout)
    case {0,1} % return fv structure for patch
        % vv = struct('Vertices',v,'Faces',f);
        vv = triangulation(f,v);
    case 2 % return vertices and faces
        vv = v; ff = f;
    otherwise
        error('Too many output variables, must be 1 or 2.');
end

end


%% Internal Functions - unit icosahedron
function [v,f] = icosahedron()
%icosahedron creates unit regular icosahedron
%   with vertices [1 0 0] and [-1 0 0].

lat=atan(1/2);
lat_norm=sqrt(1-lat^2);
% create vertices
v=[  -1    0    0;
    -lat     [-sqrt((5-sqrt(5))/8)    (1+sqrt(5))/4]*lat_norm;
    lat     [-sqrt((5+sqrt(5))/8)   (-1+sqrt(5))/4]*lat_norm;
    -lat     [-sqrt((5+sqrt(5))/8)    (1-sqrt(5))/4]*lat_norm;
    lat     [-sqrt((5-sqrt(5))/8)   -(1+sqrt(5))/4]*lat_norm;
    -lat    0    -lat_norm;
    lat      [sqrt((5-sqrt(5))/8)   -(1+sqrt(5))/4]*lat_norm;
    -lat      [sqrt((5+sqrt(5))/8)    (1-sqrt(5))/4]*lat_norm;
    lat      [sqrt((5+sqrt(5))/8)   (-1+sqrt(5))/4]*lat_norm;
    -lat      [sqrt((5-sqrt(5))/8)    (1+sqrt(5))/4]*lat_norm;
    lat    0    lat_norm;
    1    0    0;];

% create faces
f = [ 1, 2, 4;
      1, 4, 6;
      1, 6, 8;
      1, 8,10;
      1,10, 2;
     11, 2, 3;
      2, 3, 4;
      3, 4, 5;
      4, 5, 6;
      5, 6, 7;
      6, 7, 8;
      7, 8, 9;
      8, 9,10;
      9,10,11;
     10,11, 2;
     12, 3, 5;
     12, 5, 7;
     12, 7, 9;
     12, 9,11;
     12,11, 3];
end

function [v,f] = icosahedron_alt()
%icosahedron_alt creates unit regular icosahedron
%   Regular version.

t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; 
      1, t, 0; 
     -1,-t, 0; 
      1,-t, 0; 
      0,-1, t; 
      0, 1, t; 
      0,-1,-t;
      0, 1,-t; 
      t, 0,-1; 
      t, 0, 1; 
     -t, 0,-1; 
     -t, 0, 1];

% create faces
f = [ 1,12, 6; 
      1, 6, 2; 
      1, 2, 8;
      1, 8,11; 
      1,11,12;
      2, 6,10;
      6,12, 5;
     12,11, 3;
     11, 8, 7;
      8, 2, 9;
      4,10, 5;
      4, 5, 3;
      4, 3, 7;
      4, 7, 9; 
      4, 9,10;
      5,10, 6;
      3, 5,12;
      7, 3,11;
      9, 7, 8;
     10, 9, 2];

% normalise vertices to unit size
v = v./norm(v(1,:));

end

%% Internal functions - recursive subdivision
function [v,f]=subdivision(n,f,v)
    for gen = 1:n
        f_new = zeros(size(f,1)*4,3);
        for i = 1:size(f,1) % for each triangle
            tri_ind = f(i,:);
            % calculate mid points (add new points to v)
            v = [v;...
                sum(v(tri_ind([1,2]),:),1)/2;...
                sum(v(tri_ind([2,3]),:),1)/2;...
                sum(v(tri_ind([3,1]),:),1)/2];
            ind1=size(v,1)-2;
            ind2=size(v,1)-1;
            ind3=size(v,1);
            % generate new subdivided faces
            new_face = [tri_ind(1),ind1,ind3;
                tri_ind(2),ind2,ind1;
                tri_ind(3),ind3,ind2;
                ind1,ind2,ind3];
            % store new faces
            f_new(4*(i-1)+1:4*i,:) = new_face;
        end
        f = f_new; % update
    end
    
    % remove duplicate vertices
    [v,~,ind] = unique(v,'rows');
    v = v./sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);

    % reassign faces to trimmed vertex list and remove any duplicate faces
    f = unique(ind(f),'rows');
end


%% Internal Functions - frequency division inspired by Python icosphere package
function [subvertices,subfaces]=intdivision(n,f,v)
%     Subdivides mesh by adding vertices on mesh edges and faces. Each edge
%     will be divided in n segments. (For example, for n=2 one vertex is added
%     on each mesh edge, for n=3 two vertices are added on each mesh edge and
%     one vertex is added on each face.) If V and F are number of mesh vertices
%     and number of mesh faces for the input mesh, the subdivided mesh contains
%     V + F*(n+1)*(n-1)/2 vertices and F*n^2 faces.
%
%     Parameters
%     ----------
%     v : vertices list, dimension: (V,3)
%     f : faceID list,   dimension: (F,3). Index starts at 1 (MATLAB style)
%     n : subdivision frequency, integer (larger than 1 to make a change).
%
%     Returns
%     -------
%     subvertices : vertices list, dimension: (V + F*(n+1)*(n-1)/2, 3)
%     subfaces : faceID list, dimension: (F*n**2, 3)

    TR=triangulation(f,v);
    e=edges(TR);
    F=size(f,1);
    V=size(v,1);
    E=size(e,1);
    subfaces=NaN(F*n^2,3);
    subvertices=NaN(V+E*(n-1)+F*(n-1)*(n-2)/2,3);
    
    subvertices(1:size(v,1),:)=v;
    
    % Matrix for reading edgeID from verticeID
    V2E=NaN(size(v,1),size(v,1));
    for i=1:E
        V2E(e(i,1),e(i,2))=i;
        V2E(e(i,2),e(i,1))=-i;
    end
    
    template = face_template(n);
    ordering = vertex_ordering(n);
    reordered_template = ordering(template);
    
    % At this point, we have V vertices, and now we add (n-1) vertex per edge (on-edge vertices).
    w=(0:n);
    for i=1:E
        for k=2:n
            subvertices(V+(i-1)*(n-1)+k-1,:)=(w(n+2-k)*v(e(i,1),:)+w(k)*v(e(i,2),:))/n;
        end
    end
    
    % At this point we have E(n-1)+V vertices, and we add (n-1)*(n-2)/2 vertices per face (on-face vertices).
    vertices_end=E*(n-1)+V;
    r=1:n-1;
    for i=1:F
        % First, fixing connectivity. We get hold of the indices of all
        % vertices invoved in this subface: original, on-edges and on-faces.
        T = ((i-1)*(n-1)*(n-2)/2+vertices_end+1):...
            (i*(n-1)*(n-2)/2+vertices_end);
        eAB=V2E(f(i,1),f(i,2));
        eAC=V2E(f(i,1),f(i,3));
        eBC=V2E(f(i,2),f(i,3));
        AB=reverse_edge((abs(eAB)-1)*(n-1)+V+r,eAB<0);
        AC=reverse_edge((abs(eAC)-1)*(n-1)+V+r,eAC<0);
        BC=reverse_edge((abs(eBC)-1)*(n-1)+V+r,eBC<0);
        VEF=[f(i,:) AB AC BC T];
        subfaces((i-1)*n^2+1:i*n^2,:)=VEF(reordered_template);
        % Now geometry, computing positions of face vertices.
        subvertices(T,:)=inside_points(subvertices(AB,:),subvertices(AC,:));
    end
    
    subvertices=subvertices./sqrt(sum(subvertices.^2,2));
end

function vec=reverse_edge(vec,bool)
    if bool
        vec=vec(end:-1:1);
    end
end

function faces=face_template(n)
% E.g. n = 4
%              1
%             / \
%            2---3
%           / \ / \
%          4---5---6
%         / \ / \ / \
%        7---8---9---10
%       / \ / \ / \ / \
%      11--12--13--14--15

    faces=NaN(n^2,3);
    for i=0:(n-1)
        vertex0 = i*(i+1)/2;
        skip = i+1;
        for j =0:(i-1) % adding pairs of triangles, will not run for i==0
            faces(i^2+j*2+1,:)=[j+vertex0, j+vertex0+skip, j+vertex0+skip+1];
            faces(i^2+j*2+2,:)=[j+vertex0, j+vertex0+skip+1, j+vertex0+1];
        end
        % adding the last (unpaired, rightmost) triangle
        faces(skip^2,:)=[i+vertex0, i+vertex0+skip, i+vertex0+skip+1];
    end
    faces=faces+1;
end

function o=vertex_ordering(n)
% E.g. n = 4
%            1
%           / \
%          4---7
%         / \ / \
%        5---13--8
%       / \ / \ / \
%      6---14--15--9
%     / \ / \ / \ / \
%    2--10--11--12---3

    AB = 3+1:n+2;
    AC = n+3:2*n+1;
    BC = 2*n+2:3*n;
    OnFace = 3*n+1:(n+1)*(n+2)/2;
    
    o=NaN(1,(n+1)*(n+2)/2);
    o(1)=1;
    for i=1:n-1
        o(1+i*(i+1)/2)=AB(i);
        if i>=2
            o(2+i*(i+1)/2:(i+1)*(i+2)/2-1)=OnFace(((i-1)*(i-2)/2+1):(i*(i-1)/2));
        end
        o((i+1)*(i+2)/2)=AC(i);
    end
    o(n*(n+1)/2+1)=2;
    o(n*(n+1)/2+2:(n+1)*(n+2)/2-1)=BC;
    o((n+1)*(n+2)/2)=3;
end

function v=inside_points(vAB,vAC)
%     Returns coordinates of the inside                 .
%     (on-face) vertices (marked by star)              / \
%     for subdivision of the face ABC when         vAB1---vAC1
%     given coordinates of the on-edge               / \ / \
%     vertices  AB[i] and AC[i].                 vAB2---*---vAC2
%                                                  / \ / \ / \
%                                              vAB3---*---*---vAC3
%                                                / \ / \ / \ / \
%                                               .---.---.---.---.

    N=size(vAB,1);
    v=NaN(N*(N-1)/2,3);
    for i=2:N
        w=0:i;
        for k=2:i
            v((i-1)*(i-2)/2+k-1,:)=(w(i+2-k)*vAB(i,:)+w(k)*vAC(i,:))/i;
        end
    end
end