// Gmsh project created on Sat Jun  5 10:56:15 2021
SetFactory("OpenCASCADE");
//+
a = 0.2;
Point(1) = {0, 0, 0, a};
Point(2) = {1., 0, 0, a};
Point(3) = {1., 1., 0, a};
Point(4) = {0, 1., 0, a/10.};
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
//Physical Curve("BottomLine", 2) = {1};
//+
//Physical Curve("TopLine", 2) = {3};
//+
//Physical Curve("LeftLine", 2) = {4};
//+
//Physical Curve("RightLine", 2) = {2};
//+
Physical Surface("Surface", 1) = {1};
//+
Physical Point("Point", 3) = {3, 2, 1, 4};
//+
Physical Curve("Line", 2) = {4, 1, 2, 3};
