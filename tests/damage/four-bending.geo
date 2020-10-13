L  = 500.;
l  = 150.;
D  = 100;
a  = 1;  // width of the notch
d  = 30; // length of the notch

b  = 80; // widht of the refined mesh zone
x  = L-3*l;
y  = 7.5;

// mesh size, smaller finer meshes
h = 20.;
//h2 =  h/8; // mesh3
//h2 =  h/4; // mesh2
h2 =  h/2; // mesh1

Point(1) = {0,0,0,h}; 
Point(2) = {0.5*x,0,0,h};
Point(3) = {0.5*L-0.5*a,0,0,h2};
Point(4) = {0.5*L-0.5*a,d,0,h2};
Point(5) = {0.5*L+0.5*a,d,0,h2};
Point(6) = {0.5*L+0.5*a,0,0,h2};
Point(7) = {L-0.5*x,0,0,h};
Point(8) = {L,0,0,h};
Point(9) = {L,D,0,h};
Point(10)= {2*l+0.5*x,D,0,h};
Point(11)= {1*l+0.5*x,D,0,h};
Point(12)= {0.,D,0,h};
Point(13)= {0.5*L-0.5*b,0,0,h2};
Point(14)= {0.5*L+0.5*b,0,0,h2};
Point(15)= {0.5*L-0.5*b,D,0,h2};
Point(16)= {0.5*L+0.5*b,D,0,h2};
Point(100)= {0.5*L+y,0,0,h2};

Line(1) = {1,2};
Line(2) = {2,13};
Line(3) = {13,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,100};
Line(70) = {100,14};
Line(8) = {14,7};
Line(9) = {7,8};
Line(10) = {8,9};
Line(11) = {9,10};
Line(12) = {10,16};
Line(13) = {16,15};
Line(14) = {15,11};
Line(15) = {11,12};
Line(16) = {12,1};
Line(17) = {13,15};
Line(18) = {14,16};

Line Loop(1) = {1,2,17,14,15,16};
Line Loop(2) = {3,4,5,6,7,70,18,13,-17};
Line Loop(3) = {8,9,10,11,12,-18};

Plane Surface(1) = {1};
Plane Surface(2) = {2}; 
Plane Surface(3) = {3};

Physical Point(333)  = {2}; // fix point
Physical Point(444)  = {7}; // roller point
Physical Point(555)  = {10,11}; // load points
Physical Point(666)  = {100}; // monitored point

Physical Line(333)  = {16}; // left edge


Physical Surface(888) = {1,2,3}; 
