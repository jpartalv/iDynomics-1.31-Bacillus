package simulator.geometry;
import java.awt.Point;

/** 
 * @author  Jean Pierre
 */
public class CollisionEngine {

	/*public static Vector2D closestPointLineSegment( Vector2D X, Vector2D A, Vector2D B) 
	{
		Vector2D v = B.Minus( A );
		Vector2D w = X.Minus( A );
		float wDotv = w.DotProduct( v );
		float t = w.DotProduct( v ) / v.DotProduct( v );
		if (Float.isNaN(t)) t = 1;
		t = clamp( 0, 1, t );
		return A.Plus( v.Times( t ) );
	}*/

	/*public static boolean testCapsuleCapsule_old( Vector2D X, float r1, Vector2D A, Vector2D B, float r2) 
	{
		Vector2D P = closestPointLineSegment( X, A, B );
		//P.minusEquals( X ); //?
		//if( P.DotProduct( P ) <= ( r1 + r2 ) * ( r1 + r2 ) ) return P;
		if ((P.magnitude) * (P.magnitude) <=  ( r1 + r2 ) * ( r1 + r2 )) return true;
		return false;
	}*/

	public static double distance = 0;
	public static double[] intersectionPoints[] = {{0,0,0},{0,0,0}};
	/**
	 * @uml.property  name="intersectionPointsV"
	 * @uml.associationEnd  multiplicity="(0 -1)"
	 */
	public static ContinuousVector[] intersectionPointsV = new ContinuousVector[2];
	
	
	public static boolean TestCapsuleCapsule(/*Capsule capsule1, Capsule capsule2*/
			EuclideanVector capsule1, EuclideanVector capsule2, double radius1, double radius2  )
	{
		// Compute (squared) distance between the inner structures of the capsules
		//float s, t;
		//Point c1, c2;
		//float dist2 = ClosestPtSegmentSegment(capsule1.a, capsule1.b,
		//		capsule2.a, capsule2.b, s, t, c1, c2);
		double dist2 = ClosestPtSegmentSegment(capsule1, capsule2);
		distance = dist2;
			 
		// If (squared) distance smaller than (squared) sum of radii, they collide
		double radius = radius1 + radius2;//capsule1.r + capsule2.r;
		//return dist2 <= radius * radius;
		if (dist2 <= /*0*/ radius * radius) return true;
		return false;
	}

	// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
	// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
	// distance between between S1(s) and S2(t)
	public static double ClosestPtSegmentSegment(/*Point p1, Point q1, Point p2, Point q2,*/
			EuclideanVector d1, EuclideanVector d2)
	{
		
		double EPSILON = 0f; //tolerancia: might be greater than zero preferably
		double[] c1, c2;
		double s, t;
		//Vector2D d1 = q1 - p1; // Direction vector of segment S1
		//Vector2D d2 = q2 - p2; // Direction vector of segment S2
		//Vector2D r= p1-p2;
		EuclideanVector r = new EuclideanVector(d2.start,d1.start);
		double a = d1.DotProduct(d1); // Squared length of segment S1, always nonnegative
		double e = d1.DotProduct(d2); // Squared length of segment S2, always nonnegative
		double f = d2.DotProduct(r);
		// Check if either or both segments degenerate into points
		if (a <= EPSILON && e <= EPSILON) {
			// Both segments degenerate into points
			s=t= 0.0f;
			c1 = d1.start; //p1;
			c2 = d2.start; //p2;
			
			/**/
			intersectionPoints[0] = c1;
			intersectionPoints[1] = c2;
			intersectionPointsV[0] = new ContinuousVector(c1[0],c1[1],c1[2]);
			intersectionPointsV[1] = new ContinuousVector(c2[0],c2[1],c2[2]);
			/**/
			
			EuclideanVector c2c1 = new EuclideanVector(c2,c1);
			return c2c1.DotProduct(c2c1); //Dot(c1 - c2, c1 - c2);
		}
		if (a <= EPSILON) {
			// First segment degenerates into a point
			s = 0.0f;
			t=f/e; //s=0=>t= (b*s+f)/e=f/e
			t = clamp(t, 0.0f, 1.0f);
		} 
		else {
			double c = d1.DotProduct(r);
			if (e <= EPSILON) {
				// Second segment degenerates into a point
				t = 0.0f;
				s = clamp(-c / a, 0.0f, 1.0f); //t=0=>s= (b*t-c)/a=-c/a
			} 
			else {
				// The general nondegenerate case starts here
				double b = d1.DotProduct(d2);
				double denom = (a*e)-(b*b); // Always nonnegative
				// If segments not parallel, compute closest point on L1 to L2 and
				// clamp to segment S1. Else pick arbitrary s (here 0)
				if (denom != 0.0f) {
					s = clamp((b*f - c*e) / denom, 0.0f, 1.0f);
				} 
				else s = 0.0f;
				// Compute point on L2 closest to S1(s) using
				// t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
				t=(b*s+f)/e;
				// If t in [0,1] done. Else clamp t, recompute s for the new value
				// of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
				// and clamp s to [0, 1]
				if (t < 0.0f) {
					t = 0.0f;
					s = clamp(-c / a, 0.0f, 1.0f);
				} 
				else if (t > 1.0f) {
					t = 1.0f;
					s = clamp((b - c) / a, 0.0f, 1.0f);
				}
			}
		}
		c1= d1.Times(s).end; //p1+d1*s;
		c2= d2.Times(t).end; //p2+d2*t;
		
		/* extras */
		intersectionPoints[0] = c1;
		intersectionPoints[1] = c2;
		intersectionPointsV[0] = new ContinuousVector(c1[0],c1[1],c1[2]);
		intersectionPointsV[1] = new ContinuousVector(c2[0],c2[1],c2[2]);
		
		EuclideanVector c2c1 = new EuclideanVector(c2,c1);
		return c2c1.DotProduct(c2c1); //Dot(c1 - c2, c1 - c2);
	}

	/*public Vector2D testCircleCapsule( Vector2D X, float r1, Vector2D A, Vector2D B, float r2) 
	{
		Vector2D P = closestPointLineSegment( X, A, B );
		//P.minusEquals( X );
		X = P.Minus( X ); 
		if( P.DotProduct( P ) <= ( r1 + r2 ) * ( r1 + r2 ) ) return P;
		return null;
	}*/

	private static double clamp(double val, double min, double max){
		return Math.max(min, Math.min(max, val));
	}

	public static ContinuousVector[] exertForceToCapsule(ContinuousVector center, EuclideanVector length, double capsular_radius, EuclideanVector force, double mass)
	{
		//if mass is -1 take calculated volume as mass
		if (mass == -1)
		{
			mass = ((4/3)*Math.PI*Math.pow(capsular_radius, 3)) + //volume of sphere
					(Math.PI*Math.pow(capsular_radius, 2)*length.magnitude); //volume of cylinder
		}
		
		ContinuousVector[] result = new ContinuousVector[2];
		//obtain torque
		double[] _center = {center.x,center.y,center.z};
		EuclideanVector N = new  EuclideanVector(force.end,_center);
		EuclideanVector T = force.CrossProduct(N);
		
		//obtain moment of inertia
		// from: http://library.thinkquest.org/16600/advanced/rotationalinertia.shtml
		double I = (0.25*mass*Math.pow(capsular_radius, 2)) + ((mass*Math.pow(length.magnitude,2))/12);
		
		//apply rotation about torque
		EuclideanVector r = new EuclideanVector(_center, force.start);
        double netAngularVelocity = (force.CrossProduct(r).magnitude / (r.magnitude * r.magnitude));
        
        //netAngularVelocity *= I;
        
        //netAngularVelocity *= 0.1;
        
        //netAngularVelocity = 0.1;
        
        ContinuousVector[] center_head = {center,
        		new ContinuousVector(length.end[0],length.end[1],length.end[2])};
        ContinuousVector[] center_tail = {center,
        		new ContinuousVector(length.start[0],length.start[1],length.start[2])};
        ContinuousVector[] T1 = {
        		new ContinuousVector(T.start[0],T.start[1],T.start[2]),
        		new ContinuousVector(T.end[0],T.end[1],T.end[2])};

        
        
        if (T1[0].distance(T1[1]) > 0)
        {
	        result[0] = RotateVector(netAngularVelocity,center_tail,T1);
//	        if (((Double)result[0].x).equals(Double.NaN))
//	        	System.out.println("something is wrong here");
	        result[1] = RotateVector(netAngularVelocity,center_head,T1);
//			if (((Double)result[1].x).equals(Double.NaN))
//				System.out.println("something is wrong here");
        }
        
        
		return result;
	}
	
	//returns the angle and the torque
	public static double applyForceToCapsule
	(ContinuousVector center, EuclideanVector length, double capsular_radius, EuclideanVector force, double mass, EuclideanVector Torque)
	{
		//if mass is -1 take calculated volume as mass
		if (mass == -1)
		{
			mass = ((4/3)*Math.PI*Math.pow(capsular_radius, 3)) + //volume of sphere
					(Math.PI*Math.pow(capsular_radius, 2)*length.magnitude); //volume of cylinder
		}
		
		//calculate torque
		double[] _center = {center.x,center.y,center.z};
		EuclideanVector N = new  EuclideanVector(force.end,_center);
		//EuclideanVector T = force.CrossProduct(N);
		Torque = force.CrossProduct(N);
		
		//calculate moment of inertia
		// from: http://library.thinkquest.org/16600/advanced/rotationalinertia.shtml
		double I = (0.25*mass*Math.pow(capsular_radius, 2)) + ((mass*Math.pow(length.magnitude,2))/12);
		
		//calculate angle of rotation
		EuclideanVector r = new EuclideanVector(_center, force.start);
        double netAngularVelocity = (force.CrossProduct(r).magnitude / (r.magnitude * r.magnitude));
        
         
		return netAngularVelocity;
	}
	
	//return end point
	public static ContinuousVector RotateVector(double theta, 
			ContinuousVector[] v, ContinuousVector[] orientation)
    {
		//0 = start, 1 = end of vector
		//normalization of localized euclidean vector
		double norm = orientation[1].distance(orientation[0]);
        /*orientation[1].x /= norm;
        orientation[1].y /= norm;
        orientation[1].z /= norm;
		*/
		double vo_mag_x = (orientation[1].x - orientation[0].x) / norm;
        double vo_mag_y = (orientation[1].y - orientation[0].y) / norm;
        double vo_mag_z = (orientation[1].z - orientation[0].z) / norm;

		
        /*ContinuousVector vo = new ContinuousVector(orientation[1].x - orientation[0].x,
        		orientation[1].y - orientation[0].y,
        		orientation[1].z - orientation[0].z);
        */
        double mag_x = v[1].x-v[0].x;
        double mag_y = v[1].y-v[0].y;
        double mag_z = v[1].z-v[0].z;
        Quaternion Q1 = new Quaternion((double)0,mag_x,mag_y,mag_z);
        
        Quaternion Q2 = new Quaternion((float)Math.cos(theta / 2),
	            (float)(/*vo.x*/vo_mag_x * Math.sin(theta / 2)),
	            (float)(/*vo.y*/vo_mag_y * Math.sin(theta / 2)),
	            (float)(/*vo.z*/vo_mag_z * Math.sin(theta / 2)));
       Quaternion conjQ2 = Quaternion.Conjugate(Q2);

        Quaternion Q3;

        Q3 = Quaternion.Multiply(Quaternion.Multiply(Q2,Q1),conjQ2);

        //Vector result = new Vector(v.startPoint, Q3.x, Q3.y, Q3.z);
        ContinuousVector result = new ContinuousVector(v[0].x + Q3.x, v[0].y + Q3.y, v[0].z + Q3.z);
        return result;
    }
	
	//http://www.mathworks.com/matlabcentral/newsreader/view_thread/170200
	static boolean isPointInLine(ContinuousVector P1, ContinuousVector P2, ContinuousVector P )
	{
		double tol = 0.001;
		EuclideanVector A = new EuclideanVector(P,P1);
		EuclideanVector B = new EuclideanVector(P2,P1);
		EuclideanVector C = new EuclideanVector(P,P2);
		//(norm(cross(P-P1,P2-P1)) < tol) & 
		//(dot(P-P1,P2-P1) >= 0) & (dot(P-P2,P2-P1) <= 0)
		if (A.CrossProduct(B).magnitude < tol &&
			A.DotProduct(B) >= 0 && C.DotProduct(B) <= 0)
			return true;
		else 
			return false;
	}
	
	
	public static EuclideanVector dist3D_Segment_to_Segment( EuclideanVector u, EuclideanVector v)
	{
		double SMALL_NUM = 0;
	    EuclideanVector w = new EuclideanVector(v.start, u.start);
	    double    a = u.DotProduct(u);        // always >= 0
	    double    b = u.DotProduct(v);
	    double    c = v.DotProduct(v);        // always >= 0
	    double    d = u.DotProduct(w);
	    double    e = v.DotProduct(w);
	    double    D = a*c - b*b;       // always >= 0
	    double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	    double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0

	    // compute the line parameters of the two closest points
	    if (D < SMALL_NUM) { // the lines are almost parallel
	        sN = 0.0;        // force using point P0 on segment S1
	        sD = 1.0;        // to prevent possible division by 0.0 later
	        tN = e;
	        tD = c;
	    }
	    else {                // get the closest points on the infinite lines
	        sN = (b*e - c*d);
	        tN = (a*e - b*d);
	        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
	            sN = 0.0;
	            tN = e;
	            tD = c;
	        }
	        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
	            sN = sD;
	            tN = e + b;
	            tD = c;
	        }
	    }

	    if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
	        tN = 0.0;
	        // recompute sc for this edge
	        if (-d < 0.0)
	            sN = 0.0;
	        else if (-d > a)
	            sN = sD;
	        else {
	            sN = -d;
	            sD = a;
	        }
	    }
	    else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
	        tN = tD;
	        // recompute sc for this edge
	        if ((-d + b) < 0.0)
	            sN = 0;
	        else if ((-d + b) > a)
	            sN = sD;
	        else {
	            sN = (-d + b);
	            sD = a;
	        }
	    }
	    // finally do the division to get sc and tc
	    sc = (Math.abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
	    tc = (Math.abs(tN) < SMALL_NUM ? 0.0 : tN / tD);

	    // get the difference of the two closest points
	    double mag_x = w.end[0] + (u.mag_x * sc) - (tc * v.mag_x);
	    double mag_y = w.end[1] + (u.mag_y * sc) - (tc * v.mag_y);
	    double mag_z = w.end[2] + (u.mag_z * sc) - (tc * v.mag_z);
	    
	    EuclideanVector dP = new EuclideanVector(w.start, mag_x,mag_y,mag_z);  
	    
	    return dP;
	    //Vector   dP = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)

	    //return norm(dP);   // return the closest distance
	}
	
	static public boolean IntersectSegmentPlane(ContinuousVector a, ContinuousVector b, 
			/*Plane p*/
			ContinuousVector planeNormal, ContinuousVector pointInPlane, double t, ContinuousVector q)
	{
		double d = planeNormal.prodScalar(pointInPlane);
		// Compute the t value for the directed line ab intersecting the plane
		EuclideanVector ab = new EuclideanVector(a,b);// b - a;
		if (ab.magnitude == 0)
			return false;
		//t = (p.d - Dot(p.n, a)) / Dot(p.n, ab);
		//double a1 = (d - planeNormal.prodScalar(a));
		//double a2 = planeNormal.prodScalar(ab.getContinuousVector());
		
		//float     N = -dot(Pn.n, w);
		//float     D = dot(Pn.n, u);
		
		//Vector    u = S.P1 - S.P0;
	    //Vector    w = S.P0 - Pn.V0;
		
		EuclideanVector w = new EuclideanVector(pointInPlane,a);// b - a;
		
		double a1 = - planeNormal.prodScalar(w.getContinuousVector());
		double a2 = planeNormal.prodScalar(ab.getContinuousVector());
		
		t = a1 / a2 ;
		// If t in [0..1] compute and return intersection point
		if (t >= 0.0f && t <= 1.0f) {
		//q = a + t * ab;
		a.add(ab.Times(t).getContinuousVector());
		q = a;
		return true;
		}
		// Else no intersection
		return false;
	}
}
