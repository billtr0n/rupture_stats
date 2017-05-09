/* returns the dip given the fault normal vector.
	inputs
		(float) nhat1 : x1 component of fault normal vector
		(float) nhat2 : x2 component of fault normal vector
		(float) nhat3 : x3 component of fault normal vector		
	returns
		(float) dip : dip in degrees.  

	notes:  the mean fault plane is theta=90 the dip will be given
	        from [0, 180]. 
*/
float get_dip(float nhat1, float nhat2, float nhat3);

/* returns the strike given the fault normal vector.
	inputs
		(float) nhat1 : x1 component of fault normal vector
		(float) nhat3 : x3 component of fault normal vector
	returns
		(float) strike : strike w/r/t x1 axis in degrees.

	notes: simplifications are made assuming that the mean fault plane
	       normal is x3, and the surface lies in the x1, x3 plane.
*/
float get_strike(float nhat1, float nhat3);

/* return the rake 
	inputs
		(float) nhat1 : x1 component of fault normal
		(float) nhat2 : x2 component of fault normal
		(float) nhat3 : x3 component of fault normal
		(float) nhat1 : x1 component of fault normal
		(float) nhat2 : x2 component of fault normal
		(float) nhat3 : x3 component of fault normal
	returns
		(float) rake : rake

	notes: talk to steve about this before anything
*/
float get_rake(float nhat1, float nhat2, float nhat3, float su1, float su2, float su3);