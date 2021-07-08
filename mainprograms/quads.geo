//+
SetFactory("OpenCASCADE");

//+
//+
a = 0.5;
Point(1) = {0, 0, 0, a};
Point(2) = {1., 0, 0, a};
Point(3) = {1., 1., 0, a};
Point(4) = {0, 1., 0, a};
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//Transfinite Curve{2} = 4 Using Progression 2; 
//+
Transfinite Surface{1}; // deixa meus triangulos em forma de quadrilateros
//+
Recombine Surface{1};
//+
Physical Surface("plano", 1) = {1};
//+
Physical Point("fix", 3) = {1};
//+
Physical Curve("contorno", 2) = {4, 1, 2, 3};
//+
SetFactory("Built-in");
//
// o transfinite curve deixa meus quadrilateros tortos

