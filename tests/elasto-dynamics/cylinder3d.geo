//gmsh geomtery for cylinder fragmentation
R = 10.; //inner radius
H = 10.; // height
t = 1.0; // thickness

ri = R;//inner radius
ro = R+t;//outer radius

l0 = 1.0; h = l0/5;

//points
Point(1)  = {0,0,0,h};
Point(2)  = {ri,0,0,h};
Point(3)  = {0,ri,0,h};
Point(4)  = {-ri,0,0,h};
Point(5)  = {0,-ri,0,h};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(1) = {1,2,3,4};

Point(6)  = {ro,0,0,h};
Point(7)  = {0,ro,0,h};
Point(8)  = {-ro,0,0,h};
Point(9)  = {0,-ro,0,h};

Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};

Line Loop(2) = {5,6,7,8};

// -------------------
//  Surfaces
// -------------------
Plane Surface(1) = {1,2};

// -------------------
//  Volumes
// -------------------
//duplicate surfaces
Translate {0, 0, H} {
  Duplicata { Surface{1}; }
}

//extrude surfaces
Extrude {0, 0, H} {
  Surface{1};
}


// ----------------------
// Physical quantities
// ----------------------
Physical Volume("bulk") = {1}; 
Physical Surface("inner")  = {30,34,38,42};

