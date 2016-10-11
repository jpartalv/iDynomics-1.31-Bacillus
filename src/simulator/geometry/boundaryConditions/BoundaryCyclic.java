/**
 * \package simulator.geometry.boundaryConditions
 * \brief Package of boundary conditions that can be used to capture agent
 * behaviour at the boundary of the computation domain.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL 
 * "http://www.cecill.info".
 */
package simulator.geometry.boundaryConditions;

import java.util.List;
import org.jdom.Element;

import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;
import simulator.geometry.*;
import simulator.geometry.shape.*;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief BoundaryCyclic : close the system along a dimension.
 * 
 * It is computationally infeasible to simulate a micro-scale world on a 
 * macro-scale level, so oftentimes a small spatial sub-region is assumed to
 * represent the system as a whole; in this case, periodic boundaries are used
 * to remove artificial edge effects by assuming the simulated region adjoins
 * other, similar, regions. As a consequence, boundaries in some chosen
 * directions (generally for movements parallel to the substratum) are
 * periodic, which means that the solute concentrations and solute gradients
 * are constant across the boundary, and that agents traveling through one 
 * boundary will be translated to the other side of the domain
 * 
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 */
public class BoundaryCyclic extends ExternalBoundary
{
	/**
	 *  Serial version used for the serialisation of the class.
	 */
	private static final long       serialVersionUID = 1L;
	
	/**
	 * Shape object containing a construction of the boundary opposite this one.
	 */
	private IsShape _myOppShape;
	
	/**
	 * Vector that stores the intersection with the crossed boundary.
	 */
	private static ContinuousVector vectorIn;
	
	/**
	 * Used to translate a set of points to their respective points on the
	 * opposite side of the boundary.
	 */
	private static DiscreteVector translator = new DiscreteVector();
	
	/**
	 * \brief Initialises the boundary from information contained in the
	 * simulation protocol file.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param aDomain	The domain which this boundary condition is
	 * associated with.
	 * @param aBCParser	The XML tags that have declared this boundary in the
	 * protocol file.
	 */
	@Override
	public void init(Simulator aSim, Domain aDomain, XMLParser aBCParser) 
	{
		_mySide = aBCParser.getName();
		
		// in 3D, all cyclic boundaries are active
		if ( aDomain.is3D )
			activeForSolute = true;
		
		// in 2D, the x0y/xNy boundary is not active
		if( ! aDomain.is3D && _mySide.contains("x0y") )
			activeForSolute=false;
		
		readGeometry(aBCParser, aDomain);
		aDomain.addBoundary(this);
		aDomain.addBoundary(this.createOtherSide());		
	}
	
	/**
	 * \brief Read the geometry of this boundary from the protocol file and
	 * construct both this boundary and the opposite side (as this is cyclic).
	 * 
	 * @see Domain constructor.
	 * @param geometryRoot	XML tags in the protocol file that describe this
	 * boundary.
	 * @param aDomain	The domain which this boundary is associated with.
	 */
	@Override
	public void readGeometry(XMLParser geometryRoot, Domain aDomain)
	{
		List<Element> shapeList = geometryRoot.getChildrenElements("shape");
		String className;
		try
		{
			// Build first shape.
			className = "simulator.geometry.shape.";
			className += shapeList.get(0).getAttributeValue("class");
			_myShape = (IsShape) Class.forName(className).newInstance();
			_myShape.readShape(new XMLParser(shapeList.get(0)), aDomain);
			_mySide = geometryRoot.getName();
			// Build opposite side shape.
			className = "simulator.geometry.shape.";
			className += shapeList.get(1).getAttributeValue("class");
			_myOppShape = (IsShape) Class.forName(className).newInstance();
			_myOppShape.readShape(new XMLParser(shapeList.get(1)), aDomain);
		}
		catch (Exception e)
		{
			LogFile.writeError(e,
						"BoundaryCyclic.readGeometry(geometryRoot, aDomain)");
		}
	}
	
	/**
	 * \brief Method used by another which gets the indexed grid position of a
	 * continuous vector.
	 * 
	 * Some boundary conditions (e.g. BoundaryCyclic) need the input corrected
	 * due to the condition, some don't and just return the input. Maybe we'll
	 * change this at some point as to just return the input looks a bit daft
	 * - but we'll leave it here for the moment.
	 * 
	 * @param position ContinuousVector that gives the current location of an
	 * agent to check on the grid.
	 */
	@Override
	public ContinuousVector lookAt(ContinuousVector position)
	{
		// TODO Using first intersection is a quick fix.
		ContinuousVector nCC = _myShape.getIntersections(position,
									_myShape.getNormalInside()).getFirst();
		ContinuousVector bCC = getSymmetric(nCC);
		bCC.subtract(nCC);
		nCC.sendSum(bCC, position);
		return nCC;
	}
	
	/**
     * \brief Change the status of a specified LocatedGroup to note that it
     * has been identified as being outside this boundary.
     * 
     * @param aGroup LocatedGroup object which has been detected to be outside
     * the boundary.
     */
	@Override
	public void setBoundary(LocatedGroup aGroup)
	{
		aGroup.status = -1;
		// status -1 -> outside
	}
	
	/**
	 * \brief Applies the boundary condition by modifying the movement vector.
	 * 
	 * New position is orthogonal projection of the outside point on the
	 * boundary surface.
	 * 
	 * @param anAgent	The Located Agent which is attempting to cross the
	 * boundary.
	 * @param target	The target position that the agent is moving to.
	 */
	@Override
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target)
	{
		// Determine the intersection with the crossed boundary.
		// TODO Using first intersection is a quick fix.
		vectorIn = _myShape.getIntersections(anAgent.getLocation(),
											anAgent.getMovement()).getFirst();
		// Determine the remaining movement when we touch the boundary.
		target.sendDiff(target, vectorIn);
		// Apply the residual movement on the symmetric point.
		vectorIn = getSymmetric(vectorIn);
		target.add(vectorIn);
		// Compute and update the movement vector leading to this new position.
		anAgent.getMovement().sendDiff(target, anAgent.getLocation());
	}

	/**
	 * \brief Solver for the cyclic boundary condition. 
	 * 
	 * Initialises the course along the shape of the boundary, setting the
	 * values of solute near the boundary as required.
	 * 
	 * @param aSoluteGrid Grid of solute information which is to be refreshed
	 * by the solver.
	 */
	@Override
	public void refreshBoundary(SoluteGrid aSoluteGrid) 
	{
		/*
		 * For simulations in
		 * 2D: activeForSolute is false for x0y/xNy and true for x0z/xNz.
		 * 3D: activeForSolute is always true for cyclic boundaries.
		 */
		if ( activeForSolute )
		{
			// Build translator between both boundaries
			int k = (int) Math.floor(_myOppShape.getDistance(_myShape)/
												aSoluteGrid.getResolution());
			translator.set( _myOppShape.getNormalDiscrete() );
			translator.times(k - 1);
			// Initialise the course along the shape of the boundary
			_myShape.readyToFollowBoundary(aSoluteGrid);
			/*
			 * Send a point belonging to the boundary and the closest point
			 * outside the domain.
			 */
			while (_myShape.followBoundary(dcIn, dcOut, aSoluteGrid))
			{
				dcIn.add(translator);
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
			}
		}
		else
		{
			// Initialise the course along the shape of the boundary.
			_myShape.readyToFollowBoundary(aSoluteGrid);
			translator.set( _myOppShape.getNormalDiscrete() );
			translator.times(2);
			/* 
			 * Send a point belonging to the boundary and the closest point
			 * outside the domain.
			 */
			while ( _myShape.followBoundary(dcIn, dcOut, aSoluteGrid) )
			{
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
				dcOut.add(translator);
				aSoluteGrid.setValueAt(aSoluteGrid.getValueAt(dcIn), dcOut);
			}
		}
	}
	
	/**
	 * 
	 */
	@Override
	public void applyBoundary(DiscreteVector coord) 
	{
		if ( _myShape.isOutside(coord) )
		{
			DiscreteVector diff = _myShape.getRelativePosition(coord);
			coord = _myOppShape.getAbsolutePosition(diff);
		}
	}
	
	/**
	 * \brief Creates the opposite side of the cyclic boundary such that
	 * agents can 'roll around' to the other side.
	 * 
	 * @return	BoundaryCyclic object that captures the opposite side of this
	 * boundary.
	 */
	public BoundaryCyclic createOtherSide()
	{
		BoundaryCyclic out = new BoundaryCyclic();
		out.activeForSolute = this.activeForSolute;
		out._myShape = this._myOppShape;
		out._myOppShape = this._myShape;
		out._mySide = this._mySide.replaceFirst("0", "N");
		return out;
	}
	
	/**
	 * \brief Return the intersection between the opposite shape and a
	 * provided point.
	 * 
	 * @param position	A position on a boundary.
	 * @return	ContinuousVector containing the intersection between the
	 * opposite shape and provided point.
	 */
	public ContinuousVector getSymmetric(ContinuousVector position)
	{
		// TODO Using first intersection is a quick fix.
		return _myOppShape.getIntersections(position,
							_myShape.getNormalInside()).getFirst();
	}
	
	/**
	 * \brief Returns a string noting the side of the domain that this
	 * boundary condition is on.
	 * 
	 * @return String noting the side of the domain that this condition
	 * applies to (i.e. x0z, xNz, etc).
	 */
	@Override
	public String toString()
	{
		return new String("Cyclic:"+this._mySide);
	}
}