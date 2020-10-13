//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {5, 0, 0, 0.2};
//+
Point(3) = {100, 0, 0, 5};
//+
Point(4) = {100, 100, 0, 5};
//+
Point(5) = {0, 100, 0, 5};
//+
Point(6) = {0, 5, 0, 0.2};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Circle(5) = {2, 1, 6};
//+
Line Loop(6) = {4, -5, 1, 2, 3};
//+
Plane Surface(7) = {6};
//+
Physical Line("right") = {2};
//+
Physical Line("Sym_Bot") = {1};
//+
Physical Line("Sym_Left") = {4};
//
