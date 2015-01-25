module slot()
{
	render() difference()
	{
		cylinder(r=3, h=10);
		translate([0,0,3+10]) cube([3.1,3.1,20], center=true);
	}
}

module corner()
{
translate([0,0,2]) rotate([90,0,0])
render() difference()
{
	rotate([0,-10,70])
	{
		rotate([0,31.717371,18.000074]) slot();
		rotate([0,121.717476,-53.999996]) slot();
		rotate([0,121.717476,90.000000]) slot();
	}

translate([0,-12,0]) cube([20,20,20], center=true);
}
	cylinder(r=5,h=5);
}


//for(x=[0:2])
{
	//for(y=[0:1])
	{
		translate([x*18,y*18,0]) corner();
	}
}
