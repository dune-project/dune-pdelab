density = 0.5;

Point(1) = {0, 0, 0, density};
Point(2) = {15, 0, 0, density};
Point(3) = {15, 10, 0, density};
Point(4) = {0, 10, 0, density};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {7.5, 5, 0, density};
Point(6) = {10.5, 5, 0, density};
Point(7) = {4.5, 5, 0, density};
Point(8) = {7.5, 8, 0, density};
Point(9) = {7.5, 2, 0, density};
Circle(5) = {8, 5, 7};
Circle(6) = {7, 5, 9};
Circle(7) = {9, 5, 6};
Circle(8) = {6, 5, 8};
Line Loop(9) = {3, 4, 1, 2};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(11) = {9, 10};

// the whole domain:
Physical Surface(1) = {11};
// Dirichlet Boundary
Physical Line(2) = {3, 8, 5, 6, 7, 1};
// Port1
Physical Line(3) = {4};
// Port2
Physical Line(4) = {2};
