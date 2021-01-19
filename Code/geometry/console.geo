// Beam length
l=10.;
// Beam height
h=1.;
// Number of elements in x-direction
nx = 20;
// Number of elements in y-direction
ny = 10;

Point(1) = {0, 0., 0, 1.0};
Point(2) = {l, 0., 0, 1.0};
Point(3) = {l, h, 0, 1.0};
Point(4) = {0, h, 0, 1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Transfinite Line{1,3} = nx+1;
Transfinite Line{2,4} = ny+1;
Transfinite Surface{1};

Physical Line("bottom") = {1};
Physical Line("right") = {2};
Physical Line("top") = {3};
Physical Line("left") = {4};
Physical Surface("domain") = {1};
Recombine Surface{1};
