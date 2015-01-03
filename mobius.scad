module triangle(s)
{
	scale(s)
	polygon(points=[
		[-1/2, -sqrt(3)/4 + sqrt(3)/12],
		[   0, +sqrt(3)/4 + sqrt(3)/12],
		[+1/2, -sqrt(3)/4 + sqrt(3)/12],
	]);
}


thick=10;
sides=16;
radius=50;

module segment(n)
{
	translate([radius-thick/2,0,0])
	rotate([90,0,0])
	rotate([0,0,n*120/sides])
	linear_extrude(
		height=radius*3.1415*2/sides,
		twist=120/sides,
		center=true,
		slices=1
	)
	triangle(thick);
}

for(i=[1:sides])
{
	rotate([0,0,i*360/sides]) segment(i);
}

%cylinder(r=radius,height=5);

