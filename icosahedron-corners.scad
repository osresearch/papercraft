thick=7.800000;
module connector(len) {
  render() difference() {
    cylinder(r=thick/2+2, h=2*thick);
    translate([0,0,thick/2+2]) cylinder(r=thick/2, h=2*thick);
  }
  //%translate([0,0,0]) cylinder(r=thick/2, h=len*0.48);
}

module corner()
{
render() difference()
{
translate([0,0,2]) rotate([0,180,0]) {
	sphere(r=5.900000); // 11 0x7fd178c05120
	rotate([0,121.717476,-90.000000]) connector(100.000000);
	rotate([0,121.717476,-162.000000]) connector(100.000000);
	rotate([0,121.717476,54.000000]) connector(100.000000);
	rotate([0,121.717476,126.000000]) connector(100.000000);
	rotate([0,121.717476,-17.999998]) connector(100.000000);
}

translate([0,0,-5]) cube([20,20,10], center=true);
}
cylinder(r=8,h=12);
}

for(x=[0:3])
{
	for(y=[0:2])
	{
		translate([x*35,y*35,0]) corner();
	}
}
