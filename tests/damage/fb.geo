//---geometry for four bending problem ----
L  = 500.;
l  = 150.;
D  = 100;
a  = 1;  // width of the notch
d  = 10; // length of the notch

b  = 40; // width of the refined mesh zone
x  = L-3*l;
y  = 7.5;

//---- mesh size, smaller finer meshes ------------
// min dimension(Dmin) = 100 mm 
// PFM => l_0 = Dmin/100 (1mm) to Dmin/50 (2mm) =>[1.0 1.5 2.0]mm
// GED => c_0 = l_0^2/4.0 = (0.25mm^2 to 1.00mm^2) =>[0.25 0.50 1.00]mm^2
// mesh size ~ l0/5 = 0.2mm to 0.4mm
h = 20.;
hh =  0.4;

Point(1) = {0,0,0,h}; 
Point(2) = {0.5*x,0,0,h};
Point(3) = {0.5*L-0.5*a,0,0,hh};
Point(4) = {0.5*L-0.5*a,d,0,hh};
Point(5) = {0.5*L+0.5*a,d,0,hh};
Point(6) = {0.5*L+0.5*a,0,0,hh};
Point(7) = {L-0.5*x,0,0,h};
Point(8) = {L,0,0,h};
Point(9) = {L,D,0,h};
Point(10)= {2*l+0.5*x,D,0,h};
Point(11)= {1*l+0.5*x,D,0,h};
Point(12)= {0.,D,0,h};
Point(13)= {0.5*L-0.5*b,0,0,hh};
Point(14)= {0.5*L+0.5*b,0,0,hh};
Point(15)= {0.5*L-0.5*b,D,0,hh};
Point(16)= {0.5*L+0.5*b,D,0,hh};
Point(17)= {0.5*L+y,0,0,hh};

Line(1) = {1,2};
Line(2) = {2,13};
Line(3) = {13,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,17};
Line(70) = {17,14};
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

Physical Point("fix")  = {2}; // fix point
Physical Point("roller")  = {7}; // roller point
Physical Point("disp")  = {10,11}; // load points
Physical Point("monitor")  = {17}; // monitored point

Physical Surface("bulk") = {1,2,3}; 
