module tetrahedron()
{
	polyhedron(
		points=[
			[+1,0,-1/sqrt(2)],
			[-1,0,-1/sqrt(2)],
			[0,+1,+1/sqrt(2)],
			[0,-1,+1/sqrt(2)]
		],
		faces=[
			[1,0,2],
			[2,0,3],
			[0,1,3],
			[1,2,3]
		]
	);
}

scale(100) tetrahedron();
