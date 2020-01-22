--
Reference:
1. ansyshelp_11.chm
2. CALFEM - a finite element toolbox for MATLAB https://github.com/CALFEM/calfem-matlab
3. http://kis.tu.kielce.pl//mo/COLORADO_FEM/colorado/Home.html


The elements to do a piezoresistive analysis in ANSYS are: 
PLANE223, KEYOPT(1) = 101 - coupled-field 8-node quadrilateral 
SOLID226, KEYOPT(1) = 101 - coupled-field 20-node brick 
SOLID227, KEYOPT(1) = 101 - coupled-field 10-node tetrahedron 
