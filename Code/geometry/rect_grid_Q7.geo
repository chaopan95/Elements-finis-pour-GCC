// Rectangular grid of beams (dimensions L x H)
// with Nx elements in direction x and
// Ny in direction z
// inner beams are described by tag "horz" and "vert" for x and y directions
// and beams on the lateral boundaries are described by "bottom","right","top","left"

L = 1;
H = 5;
Nx = 8;
Ny = 40;
ax = L/Nx;
ay = H/Ny;

d = ax;

x = ax/2.;
y = ay/2.;
Function Cell
	p0 = newp; Point(p0) = {x, y, 0, d};
	p1 = newp; Point(p1) = {x-ax/2., y-ay/2., 0, d};
	p2 = newp; Point(p2) = {x+ax/2., y-ay/2., 0, d};

	p3 = newp; Point(p3) = {x+ax/2.,  y+ay/2.,  0.,  d} ;
	p4 = newp; Point(p4) = {x-ax/2.,  y+ay/2.,  0.,  d} ;


	l1 = newl; Line(l1) = {p1,p2};
	l2 = newl; Line(l2) = {p2,p3};
	l3 = newl; Line(l3) = {p3,p4};
	l4 = newl; Line(l4) = {p4,p1};


Return 


For j In {1:Ny}
  For i In {1:Nx}
    t = i+Nx*(j-1);
    Call Cell;
    x += ax;
    If (i==Nx)
      right[j-1]=l2;
    EndIf
    If (i==1)
      left[j-1]=l4;
    EndIf
    If (j==Ny)
      top[i-1]=l3;
    EndIf
    If (j==1)
      bot[i-1]=l1;
    EndIf
    t1 = i+Nx*(j-1); 
     If (j<Ny)
	horz[t1-1] = l3;
     EndIf
    t2 = j+Ny*(i-1); 
     If (i<Nx)
	vert[t2-1] = l2;
     EndIf
  EndFor
  x = ax/2.;
  y += ay;
EndFor


Coherence;


Physical Line("horz") = {horz[]};
Physical Line("vert") = {vert[]};
Physical Line("right") = {right[]};
Physical Line("left") = {left[]};
Physical Line("top") = {top[]};
Physical Line("bottom") = {bot[]};






















































































































