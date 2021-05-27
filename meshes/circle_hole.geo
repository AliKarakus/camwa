cl__1 = 0.25; 
cl__2 = 0.05; 
Point(1) = {-2,-2, 0, cl__1};
Point(2) = { 2,-2, 0, cl__1};
Point(3) = { 2, 2, 0, cl__1};
Point(4) = {-2, 2, 0, cl__1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(5) = {0, 0, 0, cl__2};
Point(6) = {0.5, 0, 0, cl__2};
Point(10) = {-0.5, 0, 0, cl__2};
Circle(5) = {6, 5, 10};
Circle(6) = {10, 5, 6};
Line Loop(1) = {2, 3, 4, 1};
Curve Loop(2) = {5, 6};
//+
Plane Surface(2) = {1, 2};
Physical Surface("Domain", 9) = {2};
Physical Line("Inflow", 1) = {2, 3, 4, 1,5,6};
