//create a dodecahedron by intersecting 6 boxes
module dodecahedron(height) 
{
        scale([height,height,height]) //scale by height parameter
        {
                intersection(){
                        //make a cube
                        cube([2,2,1], center = true); 
                        intersection_for(i=[0:4]) //loop i from 0 to 4, and intersect results
                        { 
                                //make a cube, rotate it 116.565 degrees around the X axis,
                                //then 72*i around the Z axis
                                rotate([0,0,72*i])
                                        rotate([116.565,0,0])
                                        cube([2,2,1], center = true); 
                        }
                }
        }
}

dodecahedron(100);
