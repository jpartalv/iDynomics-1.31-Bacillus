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

import utils.LogFile;
import utils.XMLParser;
import simulator.*;
import simulator.geometry.*;
import simulator.agent.LocatedAgent;
import simulator.agent.LocatedGroup;

/**
 * \brief Defines the bulk boundary: the concentration on the boundary is
 * fixed by a dynamic bulk, the agents crossing this line die.
 * 
 * This boundary simulates the connection to a larger bulk liquid subject to a
 * dilution process, as would be typical for a reactor. Behaviour for agents
 * does not differ from the constant concentration boundary case (and so
 * within the computational domain the boundary condition is identical), but
 * solute dynamics in the bulk compartment require an additional computational
 * step. In this step, ordinary differential equations describing the
 * reactions occurring in the biofilm and the hydraulic processes affecting 
 * the bulk liquid are solved in order to determine the bulk concentration at
 * the next time-step.
 *
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Brian Merkey (brim@env.dtu.dk, bvm@northwestern.edu), Department of
 * Engineering Sciences and Applied Mathematics, Northwestern University (USA)
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public class BoundaryBulk extends ConnectedBoundary
{
	/**
	 * Value of solute in the bulk.
	 */
	static Double	bulkValue;
	
	/**
	 * \brief Initialises the boundary from information contained in the
	 * simulation protocol file.
	 * 
	 * In this case also links the connected bulk to this boundary.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions 
	 * specified in the protocol file.
	 * @param aDomain	The domain which this boundary condition is associated
	 * with.
	 * @param aBCParser	The XML tags that have declared this boundary in the
	 * protocol file.
	 */
	@Override
	public void init(Simulator aSim, Domain aDomain, XMLParser aBCParser)
	{
		// Load the geometry of the boundary
		readGeometry(aBCParser, aDomain);		
		
		aDomain.addBoundary(this);
		
		// Load description of the connected bulk
		if ( aBCParser.isParamGiven("bulk") )
		{
			String bulkName = aBCParser.getParam("bulk");
			_connectedBulk = aSim.world.getBulk(bulkName);
		}
		else
		{
			LogFile.writeLogAlways(
					"Error! No bulk name given for BoundaryBulk\nExiting...");
			System.exit(-1);
		}
	}
	
	/**
	 * \brief Solver for the variable concentration boundary condition. 
	 * 
	 * Initialises the course along the shape of the boundary.
	 * 
	 * @param aSoluteGrid	Grid of solute information which is to be
	 * refreshed by the solver.
	 */
	@Override
	public void refreshBoundary(SoluteGrid aSoluteGrid) 
	{
		// Store the concentration in the bulk
		bulkValue = _connectedBulk.getValue(aSoluteGrid.soluteIndex);
		
		// Initialise the course along the shape of the boundary
		_myShape.readyToFollowBoundary(aSoluteGrid);

		while ( _myShape.followBoundary(dcIn, dcOut, aSoluteGrid) )
			aSoluteGrid.setValueAt(bulkValue, dcOut);
	}
	
	/**
     * \brief Change the status of a specified LocatedGroup to note that it
     * has been identified as being outside this boundary.
     * 
     * @param aGroup	LocatedGroup object which has been detected to be
     * outside the boundary.
     */
	@Override
	public void setBoundary(LocatedGroup aGroup)
	{
		aGroup.status = 3;
		// status 3 -> bulk
	}

	/**
	 * \brief Kills any agents that are crossing this boundary as they are
	 * leaving the simulated system.
	 * 
	 * @param anAgent	The LocatedAgent that has crossed the boundary.
	 * @param target	Vector of where this agent was going to be placed.
	 */
	@Override
	public void applyBoundary(LocatedAgent anAgent, ContinuousVector target) 
	{
		deadlyBoundary(anAgent, target, "overBoard");
	}


	/**
	 * \brief Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * Returns a string noting the side of the domain that this boundary condition is on
	 * 
	 * @return String noting the side of the domain that this condition applies to (i.e. x0z, xNz, etc)
	 */
	@Override
	public String toString() 
	{
		return new String("Bulk:"+this._mySide);
	}
}
