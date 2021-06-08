// Gmsh project created on Sat Jun  5 10:56:15 2021
SetFactory("OpenCASCADE");
//+
a = 0.5;
Point(1) = {0, 0, 0, a};
Point(2) = {1., 0, 0, a};
Point(3) = {1., 1., 0, a};
Point(4) = {0, 1., 0, a};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 2};
//+
Line(3) = {4, 3};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 2, -1, -4};
//+
Plane Surface(1) = {1};
Transfinite Curve{1,2,3,4} = 2 Using Progression 1;
Transfinite Surface{1};
Recombine Surface{1};
//+
Curve Loop(2) = {2, -1, -4, 3};
//+
Physical Curve("BottomLine", 5) = {1};
//+
Physical Curve("TopLine", 6) = {3};
//+
Physical Curve("LeftLine", 7) = {4};
//+
Physical Curve("RightLine", 8) = {2};
//+
Physical Surface("Surface", 1) = {1};
//+
Physical Point("Point", 11) = {3, 2, 1, 4};
//+
Physical Curve("Line", 10) = {4, 1, 2, 3};
