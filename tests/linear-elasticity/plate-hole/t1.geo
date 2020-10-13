l  = 50.;
w  = 50.;
r  = 5;


h = 2.;
h2 =  h/50;

Point(1) = {0,0,0,h}; // center of the hole

Point(2) = {r,0,0,h2};
Point(3) = {l,0,0,h};
Point(4) = {l,w,0,h};
Point(5) = {0,w,0,h};
Point(6) = {0,r,0,h2};


Circle(1) = {2,1,6};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Line Loop(1) = {-1,2,3,4,5};

Plane Surface(1) = {1};


Physical Line(333)  = {5}; // left edge
Physical Line(444)  = {2}; // bottom edge
Physical Line(555)  = {3}; // right edge

Physical Surface(888) = {1}; // substrate
