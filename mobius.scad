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
	translate([radius-thick,0,0])
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

module slice()
{
	rotate([0,0,-360/sides/2]) translate([50,-5,0]) cube([100,10,20], center=true);
	rotate([0,0,+360/sides/2]) translate([50,+5,0]) cube([100,10,20], center=true);
}

//for(i=[1:sides])
for(i=[1:sides])
{
	rotate([0,0,i*360/sides]) render() difference() 
	{
		segment(i);
		slice();
	}
}


//%cylinder(r=radius,height=5);

