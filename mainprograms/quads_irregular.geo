//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
//Transfinite Curve{1,2,3,4} = 2 Using Progression 1;
//+
Transfinite Surface{1};
//+
Recombine Surface{1};
//+
Physical Point("fix", 1) = {1};
//+
Physical Curve("contorno", 2) = {4, 1, 2, 3};
//+
Physical Surface("plano", 3) = {1};
//+
SetFactory("Built-in");
//+
Transfinite Curve {4} = 5 Using Progression 1;
//+
Transfinite Curve {3} = 5 Using Progression 1;
//+

