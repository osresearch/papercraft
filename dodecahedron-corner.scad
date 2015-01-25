module corner()
{
translate([0,0,2]) rotate([90,0,0])
render() difference()
{
	rotate([0,-10,70])
	render() difference()
	{
		sphere(r=8); // 0 0x972e1c0
		rotate([0,31.717371,18.000074]) translate([0,0,22.451374+2]) cube([3.000000,3.000000,44.902748], center=true);
		rotate([0,121.717476,-53.999996]) translate([0,0,22.451399+2]) cube([3.000000,3.000000,44.902798], center=true);
		rotate([0,121.717476,90.000000]) translate([0,0,22.451399+2]) cube([3.000000,3.000000,44.902798], center=true);
	}

translate([0,-12,0]) cube([20,20,20], center=true);
}
}


for(x=[0:2])
{
	for(y=[0:1])
	{
		translate([x*18,y*18,0]) corner();
	}
}
