thick=7.800000;
module connector(len) {
  render() difference() {
    rotate([0,0,30]) cylinder(r=thick/2+2, h=2*thick, $fa=60);
    translate([0,0,thick/2+2]) cylinder(r=thick/2, h=2*thick);
  }
  %translate([0,0,0]) cylinder(r=thick/2, h=len*0.48);
}
{
sphere(r=5.900000, $fa=60); // 1 0x7fb9f8403fb0
rotate([0,90.000000,71.999908]) connector(44.902824);
rotate([0,148.282623,-161.999924]) connector(44.902748);
rotate([0,90.000000,-35.999943]) connector(44.902859);
}
