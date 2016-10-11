package simulator.geometry.pointProcess;

import simulator.geometry.ContinuousVector;
import simulator.geometry.shape.IsShape;
import utils.ExtraMath;

/**
 * 
 * 
 *
 */
public class SweepTable
{
	/**
	 * 
	 */
	private IsShape _space;
	
	/**
	 * 
	 */
	int size;
	
	/**
	 * 
	 */
	HalfEdge[] hash;
	
	HalfEdge start, finish;
	
	/**
	 * 
	 */
	Double minValue;
	
	/**
	 * 
	 */
	Double deltaValue;
	
	private static int primary, secondary;
	
	public SweepTable(IsShape shape, int numberOfSites)
	{
		this._space = shape;
		
		primary = _space.getPrimary();
		secondary = _space.getSecondary();
		
		Double temp = 2.0 * Math.sqrt(numberOfSites + 4);
		this.size = temp.intValue();
		System.out.println("Space min "+shape.getMinPrimary()+", max "+shape.getMaxPrimary());
		this.minValue = shape.getMinPrimary();
		this.deltaValue = shape.getMaxPrimary() - minValue;
		
		this.hash = new HalfEdge[this.size];
		
		this.start = new HalfEdge(true);
		this.finish = new HalfEdge(false);
		this.start.nextNeighbor = this.finish;
		this.finish.previousNeighbor = this.start;

		this.hash[0] = start;
		this.hash[this.size - 1] = finish;
	}
	
	/**
	 * 
	 * 
	 * @param b
	 * @return
	 */
	private HalfEdge get(int b)
	{
		// If we're outside the appropriate range, return null.
		if ( b < 0 || b >= size)
			return null;
		// Find the HalfEdge corresponding to the given integer.
		HalfEdge out = hash[b];
		// If this is marked for deletion, delete it and return null.
		if ( out != null )
			if ( out.deleted )
				return (hash[b] = null);
		// Otherwise, return the HalfEdge (even if it is null).
		return out;
	}
	
	/**
	 * 
	 * 
	 * @param oldPrevious
	 * @param newNext
	 */
	public void insert(HalfEdge oldPrevious, HalfEdge newNext)
	{
		System.out.println("\nSweepTable.insert: ");
		System.out.println(newNext.toString());
		System.out.print("after");
		if ( oldPrevious.equals(this.start) )
			System.out.println(" SweepTable start");
		else
			System.out.println("\n"+oldPrevious.toString());
		
		newNext.previousNeighbor = oldPrevious;
		newNext.nextNeighbor = oldPrevious.nextNeighbor;
		(newNext.nextNeighbor).previousNeighbor = newNext;
		oldPrevious.nextNeighbor = newNext;
	}
	
	/**
	 * \brief
	 * 
	 * @param he HalfEdge to be deleted.
	 */
	public void delete(HalfEdge he)
	{
		(he.previousNeighbor).nextNeighbor = he.nextNeighbor;
		(he.nextNeighbor).previousNeighbor = he.previousNeighbor;
		he.deleted = true;
	}
	
	/**
	 * \brief find the HalfEdge immediately to the left of the given point.
	 * 
	 * @param point 
	 * @return
	 */
	public HalfEdge leftBoundary(ContinuousVector point)
	{
		/*
		 * Use hash table to get close to desired halfedge
		 */
		//System.out.println("size "+size+", value "+getValue(point)+", min "+minValue+", delta "+deltaValue);
		Double temp = size*(getValue(point)-minValue)/deltaValue;
		int bucket = temp.intValue();
		/* Ensure bucket is in the range (0, this.size - 1) */
		bucket = Math.max(Math.min(bucket, size - 1), 0);
		System.out.println("bucket "+bucket);
		HalfEdge out = get(bucket);
		/* 
		 * Starting with bucket, search backwards and forwards in the hash map
		 * to find the first non-null entry. This is our initial guess.
		 */
		if ( out == null )
			for (int i = 1; i < size; i++)
			{
				if ( (out = get(bucket - i)) != null )
					break;
				if ( (out = get(bucket + i)) != null )
					break;
			}
		/* 
		 * Linear search through the HalfEdges:
		 * 	IF: The initial guess is to the left of the Site, so keep moving
		 * 		right until the HE to the right is right of the Site.
		 *	ELSE: The initial guess is to the right of the Site, so keep
		 *		moving left until the HE is left of the Site.
		 */ 
		if (out==start || ((out!=finish) && isHErightOfPoint(out, point)))
		{
			printHalfEdge("Starting left with ", out, ", going right");
			do
			{
				out = out.nextNeighbor;
				
				printHalfEdge("\tTrying ", out, "");
				if ( out != finish )
					System.out.println("\t\tIs right of point? "+isHErightOfPoint(out, point));
			}
			while ( out != finish && isHEleftOfPoint(out, point) );
			/* This is the HalfEdge immediately the right of the point, so go
			 * left one HE.
			 */
			out = out.previousNeighbor;
			printHalfEdge("\tReturning left by one, to ", out, "");
		}
		else
		{
			printHalfEdge("Starting right with ", out, ", going left");
			// Need a do-while in case: out == rightEnd && isLeftOfSite()
			do
			{
				out = out.previousNeighbor;
				printHalfEdge("\tTrying ", out, "");
			}
			while ( out != finish && isHErightOfPoint(out, point) );
			/* Nothing more to do here - we've found the HalfEdge immediately
			 * to the left of the point.
			 */
		}
		// Update hash table.
		if (bucket > 0 && bucket < this.size - 1)
			this.hash[bucket] = out;
		return out;
	}
	
	public HalfEdge halfEdgeImmediatelyBehind(ContinuousVector point)
	{
		/*
		 * Use hash table to get close to desired halfEdge
		 */
		Double temp = size*(getValue(point)-minValue)/deltaValue;
		int bucket = temp.intValue();
		/* Ensure bucket is in the range (0, this.size - 1) */
		bucket = Math.max(Math.min(bucket, size - 1), 0);
		System.out.println("bucket "+bucket);
		HalfEdge out = get(bucket);
		/* 
		 * Starting with bucket, search backwards and forwards in the hash map
		 * to find the first non-null entry. This is our initial guess.
		 */
		if ( out == null )
			for (int i = 1; i < size; i++)
			{
				if ( (out = get(bucket - i)) != null )
					break;
				if ( (out = get(bucket + i)) != null )
					break;
			}
		/* 
		 * Linear search through the HalfEdges:
		 * 	IF: The initial guess is to the left of the Site, so keep moving
		 * 		right until the HE to the right is right of the Site.
		 *	ELSE: The initial guess is to the right of the Site, so keep
		 *		moving left until the HE is left of the Site.
		 */ 
		if ( out==start || ((out!=finish) && isPointAheadOfHE(out, point)) )
		{
			printHalfEdge("Starting behind with ", out, ", going forward");
			do
			{
				out = out.nextNeighbor;
				
				printHalfEdge("\tTrying ", out, "");
				if ( out != finish )
				{
					System.out.println("\t\tIs point ahead? "+
							isPointAheadOfHE(out, point));
				}
			}
			while ( out != finish && isPointAheadOfHE(out, point) );
			/* This is the HalfEdge immediately the right of the point, so go
			 * left one HE.
			 */
			out = out.previousNeighbor;
			printHalfEdge("\tGoing back by one, to ", out, "");
		}
		else
		{
			printHalfEdge("Starting ahead with ", out, ", going back");
			// Need a do-while in case: out == rightEnd && isLeftOfSite()
			do
			{
				out = out.previousNeighbor;
				printHalfEdge("\tTrying ", out, "");
			}
			while ( out != start && isPointBehindHE(out, point) );
			/* Nothing more to do here - we've found the HalfEdge immediately
			 * to the left of the point.
			 */
		}
		// Update hash table.
		if (bucket > 0 && bucket < this.size - 1)
			this.hash[bucket] = out;
		return out;
	}
	
	private Boolean isPointAheadOfHE(HalfEdge he, ContinuousVector point)
	{
		/* 
		 * Convert the point given, and the Site above the HalfEdge, to local
		 * coordinates.
		 */
		Double[] p = _space.convertToLocal(point);
		Double[] r = _space.convertToLocal(he.getSiteAbove());
		Double primaryDiff = p[primary] - r[primary];
		/*
		 * The easy cases: if the point is 
		 */
		Boolean rightOfSiteAbove = (primaryDiff > 0.0);
		Boolean leftOfSiteAbove = ! rightOfSiteAbove;
		
		if ( rightOfSiteAbove && he.isOutbound() )
			return true;
		if ( leftOfSiteAbove && he.isInbound() )
			return false;
		/*
		if ( point.equals(he.edge.getInnerEndPoint()) )
			return he.isInbound();
		
		if ( point.equals(he.edge.getOuterEndPoint()) )
			return he.isOutbound();
		*/
		/*
		 * Those were the simple cases: now solve the more complicated cases.
		 */
		Double secondaryDiff = p[secondary] - r[secondary];
		/*
		 * Whether or not the point given is above (true) or below (false) the
		 * point given.
		 */
		Boolean pointRightOfEdge;
		
		/*
		 * 
		 */
		ContinuousVector shadowOnEdge =
						_space.getEdgePointFromPrimary(he.edge, p[primary]);
		if ( he.edge.isNearVertical() )
		{
			/*
			 * The edge is more parallel to the secondary axis than it is to
			 * the primary axis.
			 * 
			 * In the plane, u + kv*v = K and so v = (-1/kv)*(u - K).
			 *     If kv is negative, the slope of the edge is positive and >1
			 * (-1 < kv < 0). In other words it's pointing between 12:00 and
			 * 1:30 on a clock, or North to NE on a compass.
			 *     If kv is non-negative, the slope of the edge is negative
			 * and >1 or the line is vertical (0 <= kv < 1). In other words
			 * it's pointing between 10:30 and 12:00 on a clock, or NW to
			 * North on a compass.
			 * 
			 * Unpack the edge equation to make code more readable. kv is the
			 * edge coefficient alone the secondary axis.
			 */
			Double kv = he.edge.coefficient[1];
			/*
			 * 
			 */
			Boolean fast = false;
			if ( leftOfSiteAbove && kv < 0.0 || rightOfSiteAbove && kv >= 0.0 )
			{
				/*
				 * 
				 */
				pointRightOfEdge = ( secondaryDiff >= kv * primaryDiff);
				fast = pointRightOfEdge;
			}
			else
			{
				/*
				 * 
				 */
				pointRightOfEdge = _space.compareSecondary(point, shadowOnEdge) > 0;
				if ( kv < 0.0 )
					pointRightOfEdge = ! pointRightOfEdge;
				fast = ( ! pointRightOfEdge );
			}
			/*
			 * If fast is false the problem is still not solved, and so we
			 * test further.
			 */
			if ( ! fast )
			{
				Double t1, t2;
				t1 = ExtraMath.sq(primaryDiff) - ExtraMath.sq(secondaryDiff);
				t1 *= kv;
				t2 = r[primary] - _space.getPrimary(he.getSiteBelow());
				t2 *= secondaryDiff * (1.0 + ExtraMath.sq(kv));
				t2 += 2.0 * primaryDiff * secondaryDiff;
				pointRightOfEdge = t1 < t2;
				if ( kv < 0.0 )
					pointRightOfEdge = ! pointRightOfEdge;
			}
		}
		else
		{
			/*
			 * TODO Explain why this is the way it is!
			 * 
			 * The edge is more parallel to the primary axis than it is to
			 * the secondary axis.
			 * 
			 * In the plane, ku*u + v = K and so v = -ku*u + K.
			 *     If ku is negative, the slope of the edge is positive and <1
			 * (-1 < ku < 0). In other words it's pointing between 1:30 and
			 * 3:00 on a clock, or NE to East on a compass.
			 *     If kv is non-negative, the slope of the edge is negative
			 * and <1 or the line is horizontal (0 <= ku < 1). In other words
			 * it's pointing between 3:00 and 4:30 on a clock, or East to
			 * SE on a compass.
			 */
			pointRightOfEdge = _space.distance(point, shadowOnEdge) > 
					_space.distance(he.getSiteAbove(), shadowOnEdge);
		}
		return he.isInbound() == pointRightOfEdge;
	}
	
	private Boolean isPointBehindHE(HalfEdge he, ContinuousVector point)
	{
		return ! isPointAheadOfHE(he, point);
		}
	
	/**
	 * TODO Check!
	 * 
	 * \brief 
	 * 
	 * @param halfEdge
	 * @param point
	 * @return
	 */
	private Boolean isHErightOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		return ! isHEleftOfPoint(halfEdge, point);
	}
	
	/**
	 * 
	 * @param halfEdge
	 * @param point
	 * @return
	 */
	private Boolean isHEleftOfPoint(HalfEdge halfEdge, ContinuousVector point)
	{
		/*
		 * Convert the point given, and the Site on the right of the HalfEdge,
		 * to local coordinates.
		 */
		Double[] p = _space.convertToLocal(point);
		Double[] r = _space.convertToLocal(halfEdge.getSiteAbove());
		Double primaryDiff = p[primary] - r[primary];
		/*
		 * TODO Delete above?
		 */
		Boolean rightOfSiteAbove = 
					_space.comparePrimary(point, halfEdge.getSiteAbove()) > 0;
		Boolean leftOfSiteAbove = ! rightOfSiteAbove;
		/*
		 * If the point's primary coordinate is greater than the site's, and
		 * the HE is on the left side of the Edge, then the HE is on the left
		 * of the point.
		 */
		if ( rightOfSiteAbove && halfEdge.isOutbound() )
			return true;
		/*
		 * If the point's primary coordinate is smaller than the site's, and
		 * the HE is on the right side of the Edge, then the HE is on the
		 * right of the point.
		 */
		if ( leftOfSiteAbove && halfEdge.isInbound() )
			return false;
		/*
		 * Those were the simple cases: now solve the more complicated cases.
		 */
		Double secondaryDiff = p[secondary] - r[secondary];
		Boolean pointAboveEdge = true;
		/*
		 * Unpack the edge equation to make code more readable. Names are as
		 * for edges in the plane.
		 */
		Double kv = halfEdge.edge.coefficient[1];
		/*
		 * Temporary variables to make code more readable.
		 */
		ContinuousVector pointOnEdge =
					_space.getEdgePointFromPrimary(halfEdge.edge, p[primary]);
		Double t1, t2, t3;
		if ( halfEdge.edge.isNearVertical() )
		{
			/*
			 * The edge is more parallel to the secondary axis than it is to
			 * the primary axis.
			 * 
			 * In the plane, u + kv*v = K and so v = (-1/kv)*(u - K).
			 *     If kv is negative, the slope of the edge is positive and >1
			 * (-1 < kv < 0). In other words it's pointing between 12:00 and
			 * 1:30 on a clock, or North to NW on a compass.
			 *     If kv is non-negative, the slope of the edge is negative
			 * and >1 or the line is vertical (0 <= kv < 1). In other words
			 * it's pointing between 10:30 and 12:00 on a clock, or NE to
			 * North on a compass.
			 */
			Boolean fast = false;
			if ( leftOfSiteAbove && kv < 0.0 || rightOfSiteAbove && kv >= 0.0 )
			{
				/*
				 * 
				 */
				pointAboveEdge = ( secondaryDiff >= kv * primaryDiff);
				fast = pointAboveEdge;
			}
			else
			{
				/*
				 * 
				 */
				pointAboveEdge = _space.compareSecondary(point, pointOnEdge) > 0;
				if ( kv > 0.0 )
					pointAboveEdge = ! pointAboveEdge;
				if ( ! pointAboveEdge )
					fast = true;
			}
			/*
			 * If fast is false the problem is still not solved, and so we
			 * test further.
			 */
			if ( ! fast )
			{
				t1 = ExtraMath.sq(primaryDiff) - ExtraMath.sq(secondaryDiff);
				t1 *= kv;
				t2 = r[primary] - _space.getPrimary(halfEdge.getSiteBelow());
				t2 *= secondaryDiff * (1.0 + ExtraMath.sq(kv));
				t2 += 2.0 * primaryDiff * secondaryDiff;
				pointAboveEdge = t1 < t2;
				if ( kv < 0.0 )
					pointAboveEdge = ! pointAboveEdge;
			}
		}
		else
		{
			/*
			 * TODO Explain why this is the way it is!
			 */
			pointAboveEdge = _space.distance(point, pointOnEdge) > 
					_space.distance(halfEdge.getSiteAbove(), pointOnEdge);
		}
		return halfEdge.isOutbound() == pointAboveEdge;
	}
	
	private Double getValue(ContinuousVector point)
	{
		return _space.getPrimary(point);
	}
	
	/**
	 * \brief Compare two points 
	 * 
	 * SweepTable uses the other axis to PriorityQueue.
	 * Primary-axis = x-axis in Fortune's paper.
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 * @see Voronoi.compare()
	 */
	private int compare(ContinuousVector point1, ContinuousVector point2)
	{
		Double[] p1 = _space.convertToLocal(point1);
		Double[] p2 = _space.convertToLocal(point2);
		Double temp = p1[_space.getPrimary()] - p2[_space.getPrimary()];
		int out = (int) Math.signum(temp);
		if ( out == 0 )
		{
			temp = p1[_space.getSecondary()] - p2[_space.getSecondary()];
			out = (int) Math.signum(temp);
		}
		return out;
	}
	
	public void clipAll()
	{
		// TODO
		for (HalfEdge he = start; he != finish; he = he.nextNeighbor)
			he = null;
		
	}
	
	/**
	 * Useful for testing/debugging.
	 */
	public void printSweepTable()
	{
		System.out.println("Sweep Table:");
		for ( HalfEdge he = this.start; he != null; he = he.nextNeighbor)
			System.out.println(he.toString());
	}
	
	/**
	 * Useful for testing/debugging.
	 */
	public void printHalfEdge(String preMsg, HalfEdge he, String postMsg)
	{
		String msg = preMsg;
		if ( he.equals(this.start) )
			msg += "SweepTable start";
		else if ( he.equals(this.finish) )
			msg += "SweepTable finish";
		else
			msg += he.toString();
		System.out.println(msg+postMsg);
	}
}