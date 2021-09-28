
Point(1) = {2500, 3400, 0.0, 10};
Point(2) = {0, 2900, 0.0, 10};
Point(3) = {625, 2950, 0.0, 10};
Point(4) = {1250, 3200, 0.0, 10};
Point(5) = {1600, 3150, 0.0, 10};
Point(6) = {-500, 2890, 0.0, 10};
Point(7) = {-1000, 2880, 0.0, 10};
Point(8) = {-1500, 2870, 0.0, 10};
Point(9) = {-2000, 2860, 0.0, 10};

BSpline(1) = {1, 5, 4, 3, 2, 6, 7, 8, 9};
Extrude {0, 5, 0} {
  Line{1};Layers{10};Recombine;
}
Physical Surface(6) = {5};
Physical Line(7) = {1};
Physical Line(8) = {2};
Physical Line(9) = {3};
Physical Line(10) = {4};


