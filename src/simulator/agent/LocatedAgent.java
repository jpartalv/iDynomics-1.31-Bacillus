/**
 * \package agent
 * \brief Package of utilities that create and manage agents in the simulation
 * and their participation in relevant reactions.
 * 
 * This package is part of iDynoMiCS v1.2, governed by the CeCILL license
 * under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the
 * CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 */
package simulator.agent;

import idyno.SimTimer;

import java.util.LinkedList;
import java.awt.Color;

import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;
import simulator.*;
import simulator.geometry.CollisionEngine;
import simulator.geometry.ContinuousVector;
import simulator.geometry.Domain;
import simulator.geometry.EuclideanVector;
import simulator.geometry.Quaternion;
import simulator.geometry.boundaryConditions.AllBC;

/**
 * \brief Extends ActiveAgent by adding functionality to control agent grid
 * location, agent shoving, agent death and division, and agent movement.
 *  
 * During each global timestep, agent divisions and agent growth lead to many
 * cases where neighbouring agents will overlap. A relaxation algorithm is
 * used to find iteratively the new overlap-minimising steady state
 * configuration of agent locations at the end of each timestep. 
 * 
 * @author Andreas Dötsch (andreas.doetsch@helmholtz-hzi.de), Helmholtz Centre
 * for Infection Research (Germany)
 * @author Laurent Lardon (lardonl@supagro.inra.fr), INRA, France
 * @author Sónia Martins (SCM808@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 * @author Rob Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK)
 */
public abstract class LocatedAgent extends ActiveAgent implements Cloneable 
{
	/**
	 * Temporary store of the new location this cell will move to.
	 */
	protected static ContinuousVector  _newLoc = new ContinuousVector();
	
	/**
	 * Radius of this agent.
	 */
	protected Double _radius = 0.0;
	
	/**
	 * Cell radius including any capsules.
	 */
	protected Double _totalRadius = 0.0;
	
	/**
	 * Radius at which this agent will divide.
	 */
	protected Double _myDivRadius = 0.0;
	
	/**
	 * Radius at which agent death is triggered.
	 */
	protected Double _myDeathRadius = 0.0;
	
	/**
	 * Volume of this agent.
	 */
	protected Double _volume = 0.0;
	
	/**
	 * Cell volume including any capsules.
	 */
	protected Double _totalVolume = 0.0;
	
	/**
	 * Agent position - continuous coordinates.
	 */
	protected ContinuousVector _location = new ContinuousVector();
	
	/**
	 * ContinuousVector noting the move that will be applied to the agents position.
	 */
	protected ContinuousVector _movement = new ContinuousVector();
	
	/**
	 * Direction in which this cell divides.
	 */
	protected ContinuousVector _divisionDirection = new ContinuousVector();
	
	/**
	 * List of neighbouring agents in this agent's vicinity.
	 */
	protected LinkedList<LocatedAgent> _myNeighbors = new LinkedList<LocatedAgent>();

	/**
	 * Index of the agent position on the vectorized grid.
	 */
	protected int _agentGridIndex;
	
	/**
	 * Boolean noting whether this agent is interacting with a surface (true)
	 * or not (false).
	 */
	protected Boolean _isAttached = false;

	/**
	 * Detachment priority
	 */
	public Double detPriority = 0.0;

	/**
	 * Stores the simulation time since the last division check
	 */
	public Double _timeSinceLastDivisionCheck = Double.MAX_VALUE;

	/**
	 * Distance based probability from a given neighbour (used in HGT).
	 */
	public Double _distProb = 0.0; 								
	
	/**
	 * Cumulative probability as to whether the plasmid will be transferred.
	 */
	public Double _distCumProb = 0.0; 	

	/* bacillus parameters */
	
	public double _capsular_radius = 0.389171508;//0.002948269; //aprox. 1 / 4 of the initial longitude of the bacteria	
	/**
	 * @uml.property  name="_headLocation"
	 * @uml.associationEnd  
	 */
	protected ContinuousVector         _headLocation          = new ContinuousVector();
	/**
	 * @uml.property  name="_tailLocation"
	 * @uml.associationEnd  
	 */
	protected ContinuousVector         _tailLocation          = new ContinuousVector();
	double rotationAngle = 0;
	/**
	 * @uml.property  name="torque"
	 * @uml.associationEnd  
	 */
	EuclideanVector torque = new EuclideanVector(_location,_location);


	/**
	 * \brief Constructor used to generate progenitor and initialise an object
	 * to store relevant parameters. 
	 */
	public LocatedAgent()
	{
		super();
		_speciesParam = new LocatedParam();
	}
	
	/**
	 * \brief Creates a daughter Located Agent by cloning this agent and
	 * parameter objects.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
	@Override
	@SuppressWarnings("unchecked")
	public Object clone() throws CloneNotSupportedException
	{
		LocatedAgent o = (LocatedAgent) super.clone();
		o._location = (ContinuousVector) this._location.clone();
		o._headLocation = (ContinuousVector) this._headLocation.clone();
		o._tailLocation = (ContinuousVector) this._tailLocation.clone();
		
		o._movement = (ContinuousVector) this._movement.clone();
		o._divisionDirection = (ContinuousVector)
											this._divisionDirection.clone();
		o._myNeighbors = (LinkedList<LocatedAgent>) this._myNeighbors.clone();
		o._agentGridIndex = this._agentGridIndex;
		
		return o;
	}
	
	/**
	 * \brief Create a new agent in a specified position.
	 * 
	 * @param position	Vector stating where this agent should be located.
	 */
	public void createNewAgent(ContinuousVector position) 
	{
		try 
		{
			// Get a clone of the progenitor.
			LocatedAgent baby = (LocatedAgent) sendNewAgent();
			baby.giveName();
			baby.updateSize();
			
			this._myDivRadius = getDivRadius();
			baby._myDivRadius = getDivRadius();
			baby._myDeathRadius = getDeathRadius();
			
			// Just to avoid to be in the carrier.
			// TODO Rob 13Mar2015: Is this correct?
			position.x += this._totalRadius;
			
			baby.setLocation(position);
			baby.registerBirth();
		} 
		catch (CloneNotSupportedException e) 
		{
			LogFile.writeError(e, "LocatedAgent.createNewAgent()");
		}
	}
	
	/**
	 * \brief Creates an agent of the specified species and notes the grid in
	 * which this is assigned.
	 *
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param xmlMarkUp	A species mark-up within the specified protocol file.
	 */
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser xmlMarkUp) 
	{	
		super.initFromProtocolFile(aSim, xmlMarkUp);
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
	}
	
	/**
	 * \brief Create an agent using information in a previous state or
	 * initialization file.
	 *
	 * Reads in data from the end of the singleAgentData array and then passes
	 * the remaining values onto the super class.
	 * 
	 * @param aSim	The simulation object used to simulate the conditions
	 * specified in the protocol file.
	 * @param singleAgentData	Data from the result or initialisation file
	 * that is used to recreate this agent.
	 */
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		/*
		 * Find the position to start at.
		 */
		int nValsRead = 5;
		int iDataStart = singleAgentData.length - nValsRead;
		/*
		 * This is necessary for the case when agents in a biofilm
		 * simulation are transferred into a chemostat.
		 */
		if ( Simulator.isChemostat )
			_location.reset();
		else
		{
			Double newAgentX, newAgentY, newAgentZ;
			newAgentX = Double.parseDouble(singleAgentData[iDataStart]);
			newAgentY = Double.parseDouble(singleAgentData[iDataStart+1]);
			if ( _agentGrid.is3D )
				newAgentZ = Double.parseDouble(singleAgentData[iDataStart+2]);
			else
				newAgentZ = 0.0;
			_location.set(newAgentX, newAgentY, newAgentZ);
		}
		/*
		 * Agent size.
		 */
		_radius      = Double.parseDouble(singleAgentData[iDataStart+3]);
		_totalRadius = Double.parseDouble(singleAgentData[iDataStart+4]);
		/*
		 * These are randomly generated.
		 */
		_myDivRadius = getDivRadius();
		_myDeathRadius = getDeathRadius();
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}

	
	/**
	 * \brief Called at each time step of the simulation to compute cell
	 * growth, update size, and monitor cell death and division.
	 * 
	 * Also determines whether the agent has reached the size at which it must
	 * divide.
	 */
	@Override
	protected void internalStep()
	{
		/*
		 * Compute mass growth over all compartments.
		 */
		grow();
		/*
		 * Apply this mass growth of all compounds on global radius and mass.
		 */
		updateSize();
		/*
		 * Divide if you have to.
		 */
		if ( willDivide() )
			divide();
		/*
		 * Die if you have to.
		 */
		if ( willDie() )
			die(true);
	}

	/**
	 * \brief Update the radius of the agent from the current mass (and then
	 * the volume) of the agent (EPS included).
	 */
	@Override
	public void updateSize() 
	{
		/* 
		 * Update the totalMass field (sum of the particles masses).
		 */
		updateMass();
		/*
		 * Check the mass is positive.
		 */
		if ( _totalMass < 0.0 )
			LogFile.writeLog("Warning: negative mass on agent "+sendName());
		/*
		 * Sum of (particles masses / particles density).
		 */
		updateVolume();
		/*
		 * Compute radius according to the volume.
		 */
		updateRadius();
		/*
		 * Check if by chance the agent is close enough to a support to be
		 * attached.
		 */
		if ( ! Simulator.isChemostat )
			updateAttachment();
	}

	/**
	 * \brief Captures cell division by making a clone of this agent using the
	 * makeKid method.
	 */
	public void divide()
	{
		try
		{
			makeKid();
		}
		catch (CloneNotSupportedException e)
		{
			LogFile.writeError(e, "LocatedAgent.divide()");
		}
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell division can be triggered.
	 * 
	 * @return	Boolean stating whether cell division should be triggered
	 * (true) or not (false).
	 */
	public boolean willDivide() 
	{
		/*
		 * This ensures that the checks for when to divide don't occur too
		 * often; at most they will occur at the rate of AGENTTIMESTEP.
		 */
		_timeSinceLastDivisionCheck += SimTimer.getCurrentTimeStep();
		if ( _timeSinceLastDivisionCheck < _agentGrid.getAgentTimeStep() )
			return false;
		_timeSinceLastDivisionCheck = 0.0;
		/*
		 * At this point we will actually check whether to divide.
		 */
		return getRadius(false) > this._myDivRadius;
	}

	/**
	 * \brief Determines whether or not a cell has reached the radius where
	 * cell death can be triggered.
	 * 
	 * @return	Boolean stating whether cell death should be triggered (true)
	 * or not (false).
	 */
	public boolean willDie()
	{
		return (_totalMass < 0.0) || (getRadius(false) <= _myDeathRadius);
	}
	
	/**
	 * \brief With it determined that cell division will occur, create a new
	 * agent from the existing one.
	 * 
	 * @throws CloneNotSupportedException Thrown if the agent cannot be cloned.
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		/*
		 * Create the new instance.
		 */
		LocatedAgent baby = (LocatedAgent) sendNewAgent();
		/*
		 * These are all generated randomly.
		 */
		this._myDivRadius = getDivRadius();
		baby._myDivRadius = getDivRadius();
		baby._myDeathRadius = getDeathRadius();
		/*
		 * Update the lineage.
		 */
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * In a chemostat, the daughter cells remain with the coordinates of
		 * their progenitor. Otherwise, compute movement to apply to both
		 * cells and apply it.
		 */
		/*if ( ! Simulator.isChemostat )
		{
			setDivisionDirection(getInteractDistance(baby)/2);
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}*/
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		baby._netVolumeRate = 0.0;
	}

	/**
	 * \brief On agent division, divides the mass between the old and new
	 * agent, at a specified fraction.
	 * 
	 * @param baby	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The fraction of this agents mass that should be
	 * transferred to the new agent.
	 */
	public void divideCompounds(LocatedAgent baby, Double babyMassFrac)
	{
		
		this._radius *= babyMassFrac;
		baby._radius *= 1- babyMassFrac;
		
		//relocate centers
		EuclideanVector orientation = new EuclideanVector(_tailLocation,_headLocation);
		orientation = orientation.Normalize();
		orientation = orientation.Times(_radius /*+ _capsular_radius*/);
		
		_location.add(orientation.getContinuousVector());
		
		orientation = new EuclideanVector(baby._tailLocation,baby._headLocation);
		orientation = orientation.Normalize();
		orientation = orientation.Times(baby._radius /*+ _capsular_radius*/);
		baby._location.subtract(orientation.getContinuousVector());
		
		//add slight Variation on orientation
		//ContinuousVector end = new ContinuousVector(ExtraMath.random.nextFloat(),ExtraMath.random.nextFloat(),ExtraMath.random.nextFloat());
		
		double proportion = 0.2;
		
		if (this._agentGrid.is3D)
		{
			torque.mag_y = ExtraMath.random.nextDouble();
			torque.mag_z = ExtraMath.random.nextDouble();
		}	
		else
		{
			if (ExtraMath.random.nextDouble() > 0.5)
				torque.mag_z = 1;
			else
				torque.mag_z = -1;
		}
		
		rotationAngle = ExtraMath.random.nextDouble() * proportion; //multiply by a PROTOCOL magnitude
		
		
		/*
		 * Apply babyMassFrac.
		 */
		for (int i = 0; i<particleMass.length; i++)
		{
			baby.particleMass[i] *= babyMassFrac;
			this.particleMass[i] *= 1-babyMassFrac;
		}
		/*
		 * Update radius, mass, volumes and growth rates.
		 */
		updateSize();
		baby.updateSize();
		updateGrowthRates();
		baby.updateGrowthRates();
	}

	/**
	 * \brief On agent division, transfers biomass and EPS between the old and
	 * new agent, at a specified ratio.
	 * 
	 * @param baby	The new agent, which is inheriting mass.
	 * @param babyMassFrac	The ratio of the biomass/EPS that should be 
	 * transferred to the new agent.
	 */
	public void transferCompounds(LocatedAgent baby, Double babyMassFrac)
	{
		Double massToTransfer;
		for (int i = 0; i<particleMass.length; i++)
		{
			massToTransfer = this.particleMass[i] * babyMassFrac;
			baby.particleMass[i] += massToTransfer;
			this.particleMass[i] -= massToTransfer;
		}
		/*
		 * Update radius, mass and volumes.
		 */
		updateSize();
		baby.updateSize();
	}
	
	/**
	 * \brief Set the movement vector that states where to put a newly-created
	 * particle.
	 * 
	 * @param distance	Distance between the this agent and the new agent.
	 */
	public void setDivisionDirection(Double distance)
	{
		Double phi, theta;
		phi = 2*Math.PI*ExtraMath.getUniRandAngle();//ExtraMath.getUniRandAngle();
		theta = 2*Math.PI*ExtraMath.getUniRandAngle();//ExtraMath.getUniRandAngle();
		_divisionDirection.x = distance * Math.sin(phi) * Math.cos(theta);
		_divisionDirection.y = distance * Math.sin(phi) * Math.sin(theta);
		if ( _agentGrid.is3D )
			_divisionDirection.z = distance * Math.cos(phi);
		else
			_divisionDirection.z = 0.0;
	}

	/* ______________________ SHOVING ___________________________________ */

	/**
	 * \brief Models a mechanical interaction between two located agents.
	 * 
	 * Implemented by extending classes (LocatedAgent).
	 * 
	 * @param MUTUAL	Whether movement is shared between two agents or
	 * applied only to this one. Set in the protocol file.
	 * @return	The move to be applied once the shoving or pull calculations
	 * have been performed.
	 */
	@Override
	public Double interact(boolean MUTUAL)
	{
		move();
		
		/*
		 * Rebuild your neighbourhood.
		 */
		getPotentialShovers(getInteractDistance());
		for ( LocatedAgent neighbour : _myNeighbors )
			addPushMovement(neighbour, MUTUAL);
		_myNeighbors.clear();
		return move();
	}

	/**
	 * \brief Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector.
	 * 
	 * Mutual shoving : The movement by shoving of an agent is calculated based on the cell overlap and added to the agents movement vector. 
	 * Both agents are moved of half the overlapping distance in opposite directions.
	 * 
	 * @param aNeighbour	 Reference to the potentially shoving neighbour
	 * @param isMutual	Whether movement is shared between two agents or applied only to this one
	 * @param gain	Double noting change in position
	 * @return Boolean stating whether shoving is detected (true) or not (false)
	 */
	public void addPushMovement(LocatedAgent aNeighbor, boolean isMutual)
	{
		/*
		 * Cannot push oneself!
		 */
		if ( aNeighbor == this )
			return;

		/* verify intersection of capsules */
		EuclideanVector bactMe = new EuclideanVector(_tailLocation,_headLocation);
		EuclideanVector bactHim = new EuclideanVector(aNeighbor._tailLocation,aNeighbor._headLocation);
		
		double dotProduct = bactMe.DotProduct(bactHim);
   		double angle = dotProduct / (bactMe.magnitude * bactHim.magnitude); 
    	
   		// if vector are in opposite direction the algorithm fails.
   		if (angle < 0)
   			bactHim = new EuclideanVector(aNeighbor._headLocation,aNeighbor._tailLocation);	
		
   		boolean theyIntersect = CollisionEngine.TestCapsuleCapsule(bactMe, bactHim, _capsular_radius, aNeighbor._capsular_radius);
		if (theyIntersect)
		{
			//calculate translation and rotation
			EuclideanVector forceMe = new EuclideanVector(CollisionEngine.intersectionPointsV[0],CollisionEngine.intersectionPointsV[1]);
			double newMag = forceMe.magnitude - (2* _capsular_radius);
			forceMe = forceMe.Normalize();
			forceMe = new EuclideanVector(forceMe.start,forceMe.mag_x * newMag, 
			forceMe.mag_y * newMag, forceMe.mag_z * newMag);
			
			//forceMe.Times(0.5f);
				
			double[] _center = {_location.x,_location.y,_location.z};
			EuclideanVector N = new  EuclideanVector(forceMe.end,_center);
			EuclideanVector T = forceMe.CrossProduct(N);;
			this.rotationAngle += CollisionEngine.applyForceToCapsule(
					this._location, new EuclideanVector(_tailLocation,_headLocation),
					_capsular_radius, forceMe, -1, null);
			torque = torque.Plus(T);
			//System.out.println(this.rotationAngle+"");
						
			if (isMutual) {
				forceMe.Times(0.5f);
				this.rotationAngle *= 0.5;
				this._movement.add(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);

				aNeighbor.rotationAngle -= rotationAngle; 
				aNeighbor.torque = /*T;*/aNeighbor.torque.Minus(T);

				aNeighbor.rotationAngle *= 0.5;
				aNeighbor._movement.subtract(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);
			} else {
				this._movement.add(forceMe.mag_x,forceMe.mag_y,forceMe.mag_z);
			}
		}
		
		return;
	}

	/**
	 * \brief Computes the shortest vector between this agent and a position
	 * given as a ContinuousVector. Assumes cyclic boundaries.
	 * 
	 * If the vector is all zero's, returns a vector of random direction and
	 * length = 0.01 * radius.
	 * 
	 * TODO Can we do this without assuming cyclic boundaries? I.e. actually
	 * check..
	 * 
	 * @param position	ContinuousVector of position to calculate distance to.
	 * @return The shortest movement vector to go from a to b, taking into
	 * account the cyclic boundary.
	 * @see addOverlapMovement
	 * @see addPullMovement works in 2 and 3D
	 */
	public ContinuousVector computeDifferenceVector(ContinuousVector position)
	{
		Double gridLength;
		ContinuousVector diff = new ContinuousVector();
		diff.sendDiff(_location, position);
		/*
		 * Check periodicity in X.
		 */
		gridLength = _species.domain.length_X;
		if ( Math.abs(diff.x) > 0.5 * gridLength )
			diff.x -= Math.signum(diff.x) * gridLength;
		/*
		 * Check periodicity in Y.
		 */
		gridLength = _species.domain.length_Y;
		if ( Math.abs(diff.y) > 0.5 * gridLength )
			diff.y -= Math.signum(diff.y) * gridLength;
		/*
		 * Check periodicity in Z.
		 */
		if (_agentGrid.is3D)
		{
			gridLength = _species.domain.length_Z;
			if (Math.abs(diff.z) > 0.5 * gridLength)
				diff.z -= Math.signum(diff.z) * gridLength;
		}
		/*
		 * If this is a zero vector, give it random direction and a norm of
		 * 0.01 * radius.
		 */
		if ( diff.isZero() )
		{
			diff.alea(_agentGrid.is3D);
			diff.normalizeVector(0.01*_radius);
		}
		return diff;
	}
	
	/**
	 * 
	 * @param aLoc
	 * @return
	 */
	public ContinuousVector computeDifferenceVector(LocatedAgent aLoc)
	{
		return computeDifferenceVector(aLoc._location);
	}
	
	/**
	 * \brief Find neighbouring agents in a range around you.
	 * 
	 * @param radius	The distance to search around the agent location.
	 */
	public void getPotentialShovers(Double radius)
	{
		_agentGrid.getPotentialShovers(_agentGridIndex, radius, _myNeighbors);
	}

	/**
	 * \brief Pick a random neighbour from the _myNeigbors collection.
	 * 
	 * If used multiple times without changing _myNeighbours, this will be
	 * random selection WITH repetition.
	 * 
	 * @return	A randomly picked neighbour (LocatedAgent object) from the
	 * list of neighbours.
	 */
	public LocatedAgent pickNeighbor()
	{
		if ( _myNeighbors.isEmpty() )
			return null;
		return _myNeighbors.get(ExtraMath.getUniRandInt(_myNeighbors.size()));
	}

	/**
	 * \brief Find siblings of this agent in the immediate surroundings.
	 * 
	 * @param indexSpecies	The index used to reference this species in the
	 * simulation dictionary.
	 */
	public void findCloseSiblings(int indexSpecies) 
	{
		Double shoveDist;
		LocatedAgent aNb;
		/*
		 * Find and count neighbours.
		 */
		getPotentialShovers(getInteractDistance());
		int nNb = _myNeighbors.size();
		/*
		 * Loop through them, only re-appending them to the neighbour list
		 * if they are: (1) different to this agent, (2) the same species as 
		 * this agent, and (3) close enough to this agent.  
		 */
		for ( int iNb = 0; iNb < nNb; iNb++ )
		{
			aNb = _myNeighbors.removeFirst();
			if ( aNb == this || indexSpecies != aNb.speciesIndex)
				continue;
			shoveDist = 2 * (getShoveRadius() + aNb.getShoveRadius());
			if ( getDistance(aNb) <= shoveDist )
				_myNeighbors.addLast(aNb);
		}
	}

	/**
	 * \brief With the agent move calculated, apply this movement, taking care
	 * to respect boundary conditions.
	 * 
	 * @return Distance moved relative to total radius.
	 */
	@Override
	public Double move()
	{
		/*
		 * Check the movement is valid.
		 */
		if ( ! _movement.isValid() )
		{
			LogFile.writeLog("Incorrect movement coordinates");
			_movement.reset();
		}
		/*
		 * Check we're not trying to move in the Z direction in 2D.
		 */
		if ( !(_agentGrid.is3D) && !(_movement.z.equals(0.0)) )
		{
			_movement.z = 0.0;
			_movement.reset();
			LogFile.writeLog("Agent tried to move in Z direction!");
		}
		
		//even if there is no movement there could be intersection with the boundaries
		//from growth
		checkBoundariesTailHead();
		_headLocation.add(_movement);
		_tailLocation.add(_movement);
		_location.add(_movement);
		
		//sometimes the capsule is still in a non valid position
		
		AllBC aBoundary = getDomain().testCrossedBoundary(_location);
		boolean testCenter = (aBoundary!=null);
		
		if (testCenter)
		{
			EuclideanVector force = 
			 new EuclideanVector(_location,aBoundary.getOrthoProj(_location));

			
			 //we need a vector pointing inside with the size of the capsular radius
			 ContinuousVector forceNormal = new ContinuousVector(0d, 0d, 0d); 
			 forceNormal = aBoundary.getShape().getNormalInside();
			 forceNormal.normalizeVector();
			 forceNormal.times(_capsular_radius);
			 force.mag_x += forceNormal.x;
			 force.mag_y += forceNormal.y;
			 force.mag_z += forceNormal.z;
			 
			 _headLocation.add(force.getContinuousVector());
			 _tailLocation.add(force.getContinuousVector());
			 _location.add(force.getContinuousVector());
			 _movement.add(force.getContinuousVector());
		}

		/*
		 * Now apply the movement.
		 */
		_agentGrid.registerMove(this);
		/*
		 * Calculate how far we've traveled relative to the total radius.
		 */
		Double delta = _movement.norm();
		_movement.reset();
		return delta/_totalRadius;
	}

	/**
	 * \brief Used by the move method to determine if an agent's move crosses
	 * any of the domain's boundaries.
	 */
	public ContinuousVector getVerifiedLocationFromMovement(ContinuousVector movement) {
		// Search a boundary which will be crossed
		ContinuousVector newLoc = new ContinuousVector(_location);
		newLoc.add(movement);
		return getVerifiedLocation(newLoc);
	}
	
	public ContinuousVector getVerifiedLocation(ContinuousVector location) {
		AllBC aBoundary = getDomain().testCrossedBoundary(getRadius(true),location);
		int nDim = (_agentGrid.is3D ? 3 : 2);
		int counter = 0;

		/*
		 * Test all boundaries and apply corrections according to crossed
		 * boundaries.
		 */
		while (aBoundary != null || (counter > nDim))
		{
			counter++;
			aBoundary.applyBoundary(this, location);
			aBoundary = getDomain().testCrossedBoundary(getRadius(true),location);
		}
		return location;
	}

	
	/**
	 * \brief Add the reacting concentration of an agent to the received grid
	 * 
	 * Add the reacting concentration of an agent to the received grid
	 * 
	 * @param aSpG	Spatial grid used to sum catalysing mass
	 * @param catalystIndex	Index of the compartment of the cell supporting the reaction
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG, int catalystIndex)
	{
		if (isDead)
			return;

		Double value = particleMass[catalystIndex]/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total concentration of an agent on received grid
	 * 
	 * Add the total concentration of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum total mass
	 */
	@Override
	public void fitMassOnGrid(SpatialGrid aSpG) 
	{
		if (isDead)
			return;

		Double value = _totalMass/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		aSpG.addValueAt(value, _location);
	}

	/**
	 * \brief Add the total volume rate of an agent on received grid
	 * 
	 * Add the total volume rate of an agent on received grid
	 * 
	 * @param aSpG	Spatial grid used to sum volume
	 */
	public void fitVolRateOnGrid(SpatialGrid aSpG)
	{
		Double value = _netVolumeRate/aSpG.getVoxelVolume();
		if ( ! Double.isFinite(value) )
			value = 0.0;
		try
		{
			aSpG.addValueAt(value, _location);
		}
		catch (ArrayIndexOutOfBoundsException e)
		{
			LogFile.writeLogAlways("Could not put LocatedAgent mass on grid");
			LogFile.writeLogAlways("Problem with location "
													+_location.toString());
			System.exit(-1);
		}
	}

	/**
	 * \brief Add the reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * Add the total reaction/growth rate of an agent on received grid, for a specified reaction
	 * 
	 * @param aRateGrid	Spatial grid used to store total reaction rate
	 * @param reactionIndex	Index of this declared reaction in the simulation dictionary
	 */
	@Override
	public void fitReacRateOnGrid(SpatialGrid aRateGrid, int reactionIndex)
	{
		if (isDead)
			return;
		
		// growthRate is in [fgX.hr-1] so convert to concentration:
		// [fgX.um-3.hr-1 = gX.L-1.hr-1]
		Double value = growthRate[reactionIndex]/aRateGrid.getVoxelVolume();

		if ( ! Double.isFinite(value) )
			value = 0.0;

		aRateGrid.addValueAt(value, _location);
	}

	/* _______________ FILE OUTPUT _____________________ */

	/**
	 * \brief Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * Used in creation of results files - specifies the header of the columns of output information for this agent
	 * 
	 * @return	String specifying the header of each column of results associated with this agent
	 */
	@Override
	public StringBuffer sendHeader()
	{
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = super.sendHeader();
		
		// location info and radius
		tempString.append(",locationX,locationY,locationZ,radius,totalRadius");
		
		return tempString;
	}

	/**
	 * \brief Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * Used in creation of results files - creates an output string of information generated on this particular agent
	 * 
	 * @return	String containing results associated with this agent
	 */
	@Override
	public StringBuffer writeOutput()
	{
		// write the data matching the header file
		StringBuffer tempString = super.writeOutput();
		
		// location info and radius
		tempString.append(","+_location.x+","+_location.y+","+_location.z+",");
		tempString.append(_radius+","+_totalRadius);
		
		return tempString;
	}

	/* _______________ RADIUS, MASS AND VOLUME _____________________ */

	/**
	 * \brief Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 * 
	 * Compute the volume on the basis of the mass and density of different compounds defined in the cell
	 */
	public void updateVolume()
	{
		_volume = 0.0;
		for (int i = 0; i<particleMass.length; i++) {
			_volume += particleMass[i]/getSpeciesParam().particleDensity[i];
		}
		_totalVolume = _volume;
	}

	/**
	 * \brief Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 * 
	 * Compute the radius on the basis of the volume The radius evolution is stored in deltaRadius (used for shrinking)
	 */
	public void updateRadius() {

		//sonia:chemostat 22.02.2010
		if (_radius == 0)
		{
			if(Simulator.isChemostat || _species.domain.is3D){
				_radius = ExtraMath.radiusOfASphere(_volume);
				_totalRadius = ExtraMath.radiusOfASphere(_totalVolume);
			}else{
				_radius = ExtraMath.radiusOfACylinder(_volume,
						_species.domain.length_Z);
				_totalRadius = ExtraMath.radiusOfACylinder(_totalVolume,
						_species.domain.length_Z);
			}
		}
		else
		{
			double reaction = ((LocatedParam) _speciesParam).reactionKinetic[0][0] ;
			_radius *= 1 + ((1 - reaction) / 60) ;
		}
		//also update length 
		if (_radius > 0)
		{
			
			//find head and tail locations given the current radius
			//first get unit vector from current head and tail then multiply by radius
			double magnitude = _headLocation.distance(_location);
			if (magnitude == 0)
			{
				randomizeOrientation();
				 magnitude = _headLocation.distance(_location);
			}
			
			double magX = ((_headLocation.x - _location.x) / magnitude) * (_radius - _capsular_radius) ;
			double magY = ((_headLocation.y - _location.y) / magnitude) * (_radius - _capsular_radius) ;
			double magZ = ((_headLocation.z - _location.z) / magnitude) * (_radius - _capsular_radius) ;
			_headLocation.x = _location.x + magX;
			_headLocation.y = _location.y + magY;
			_headLocation.z = _location.z + magZ;
			_tailLocation.x = _location.x - magX;
			_tailLocation.y = _location.y - magY;
			_tailLocation.z = _location.z - magZ;
		}
	}

	/**
	 * \brief Update the attachment, checking if this agent is close enough to
	 * any boundaries.
	 * 
	 * TODO Rob 13Mar2015: Where does this 3 come from?!
	 * 
	 * @return	Boundary that has been crossed.
	 */
	public void updateAttachment()
	{
		for (AllBC aBoundary : getDomain().getAllBoundaries())
			if ( aBoundary.isSupport() &&
						aBoundary.getDistance(_location) <= 3 * _totalRadius )
			{
				_isAttached = true;
				return;
			}
	}
	
	/**
	 * \brief Add movement to the ContinuousVector storing the agents move.
	 * 
	 * @param aMove	ContinuousVector to add to the movement vector.
	 */
	public void addMovement(ContinuousVector aMove)
	{
		this._movement.add(aMove);
	}
	
	/**
	 * \brief Return the set of parameters associated with this agent
	 * (LocatedParam object).
	 * 
	 * @return LocatedParam object of parameters associated with this agent.
	 */
	@Override
	public LocatedParam getSpeciesParam()
	{
		return (LocatedParam) _speciesParam;
	}
	
	/**
	 * \brief Return the volume of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the volume of this agent.
	 */
	public Double getVolume(boolean withCapsule)
	{
		return withCapsule ? _totalVolume : _volume;
	}
	
	/**
	 * \brief Return the radius of this agent, with or without the capsule.
	 * 
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the radius of this agent
	 */
	public Double getRadius(boolean withCapsule)
	{
		return (withCapsule ? _totalRadius : _radius);
	}
	
	/**
	 * \brief Return the mass of this agent, with or without the capsule.
	 *  
	 * @param withCapsule	Boolean noting whether any capsule should be
	 * included in this calculation.
	 * @return	Double specifying the mass of this agent.
	 */
	public Double getMass(Boolean withCapsule)
	{
		return (withCapsule ? _totalMass : _totalMass);
	}
	
	/**
	 * \brief Report whether this cell has any EPS.
	 * 
	 * @return	Boolean noting whether this cell has any EPS.
	 */
	public Boolean hasEPS() 
	{
		return false;
	}
	
	/**
	 * \brief Return the shove factor to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove factor that will be applied.
	 */
	public Double getShoveFactor()
	{
		return ((LocatedParam) _speciesParam).shoveFactor;
	}

	/**
	 * \brief Return the shove radius to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove radius that will be applied.
	 */
	public Double getShoveRadius()
	{
		return _totalRadius * getShoveFactor();
	}
	
	/**
	 * \brief Return the shove limit to be used in shoving for this
	 * species of agent.
	 * 
	 * @return	Double specifying the shove limit that will be applied.
	 */
	public Double getShoveLimit()
	{
		return ((LocatedParam) _speciesParam).shoveLimit;
	}
	
	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * for this species of agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance()
	{
		return getInteractDistance(this);
	}
	
	/**
	 * \brief Return the shoving interaction distance to be used in shoving
	 * against a specified agent.
	 * 
	 * @return	Double specifying the shoving interaction distance that will
	 * be applied.
	 */
	public Double getInteractDistance(LocatedAgent aLoc)
	{
		return getShoveRadius() + aLoc.getShoveRadius() + getShoveLimit();
	}
	
	/**
	 * \brief Return the fraction of mass that is transferred to the new agent
	 * on cell division.
	 * 
	 * @return	Double stating the fraction of mass that is transferred to the
	 * new agent on cell division.
	 */
	public Double getBabyMassFrac()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().babyMassFrac,
											getSpeciesParam().babyMassFracCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell division is triggered.
	 * 
	 * @return	Double stating the agent radius at which cell division is
	 * triggered.
	 */
	public Double getDivRadius()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().divRadius,
											getSpeciesParam().divRadiusCV);
	}
	
	/**
	 * \brief Return the agent radius at which cell death is triggered
	 * 
	 * @return	Double stating the agent radius at which cell death is triggered
	 */
	public Double getDeathRadius()
	{
		return ExtraMath.deviateFromCV(getSpeciesParam().deathRadius,
											getSpeciesParam().deathRadiusCV);
	}
	
	/**
	 * \brief Report if this agent is attached to a surface.
	 * 
	 * @return Boolean noting whether the agent is attached to a surface.
	 */
	public Boolean isAttached()
	{
		return _isAttached;
	}
	
	/**
	 * \brief Return the active fraction of this agent.
	 * 
	 * @return	Double value stating the active fraction of this agent.
	 */
	public Double getActiveFrac()
	{
		return 1.0;
	}
	
	/**
	 * \brief Return the color assigned to this agent in POV-Ray output.
	 * 
	 * TODO Rob 13Mar2015: Consider deleting as part of move away from POV-Ray.
	 * 
	 * @return	Colour assigned to this agent as specified in the protocol
	 * file.
	 */
	public Color getColor()
	{
		return _species.color;
	}

	/**
	 * \brief Return the colour assigned to any capsules contained in this
	 * agent in POV-Ray output.
	 * 
	 * @return	Colour assigned to this agent capsules as specified in the
	 * protocol file.
	 */
	public Color getColorCapsule()
	{
		return Color.green;
	}

	/**
	 * \brief Return the location of this agent.
	 * 
	 * @return	ContinuousVector stating the location of this agent.
	 */
	public ContinuousVector getLocation()
	{
		return _location;
	}

	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * Comparator used by AgentContainer.erodeBorder()
	 * @author Rob Clegg
	 */
	public static class detPriorityComparator implements java.util.Comparator<Object>
	{
		@Override
		public int compare(Object b1, Object b2)
		{
			Double f1 = ((LocatedAgent) b1).detPriority;
			Double f2 = ((LocatedAgent) b2).detPriority;
			return (int) Math.signum(f1 - f2);
		}
	}

	/**
	 * \brief Comparator used by AgentContainer.erodeBorder()
	 * 
	 * @author Rob Clegg
	 */
	public static class totalMassComparator implements java.util.Comparator<Object>
	{
		@Override
		public int compare(Object b1, Object b2)
		{
			Double f1 = ((LocatedAgent) b1)._totalMass;
			Double f2 = ((LocatedAgent) b2)._totalMass;
			return (int) Math.signum(f1 - f2);
		}
	}
	
	/**
	 * \brief Return the distance from this agent to a ContinuousVector.
	 * 
	 * @param position	ContinuousVector to find distance to.
	 * @return distance between this agent and cV (assuming cyclic boundaries).
	 */
	public Double getDistance(ContinuousVector position)
	{
		return computeDifferenceVector(position).norm();
	}
	
	/**
	 * \brief Return the distance from this agent to another.
	 * 
	 * @param aLoc	LocatedAgent to find distance to.
	 * @return Distance from this agent to that given (assuming cyclic
	 * boundaries).
	 */
	public Double getDistance(LocatedAgent aLoc)
	{
		return getDistance(aLoc._location);
	}

	/**
	 * \brief Set the location of this agent to the supplied continuous vector.
	 * 
	 * @param cc	Location which this agent should be assigned to.
	 */
	public void setLocation(ContinuousVector cc) 
	{
		// In a chemostat set the location of the newborns to zero.
		if ( Simulator.isChemostat )
			_location.reset();
		else
		{
			_headLocation.x += (cc.x - _location.x);
			_headLocation.y += (cc.y - _location.y);
			_headLocation.z += (cc.z - _location.z);
			_tailLocation.x += (cc.x - _location.x);
			_tailLocation.y += (cc.y - _location.y);
			_tailLocation.z += (cc.z - _location.z);
			
			_location.set(cc);
			
		}
	}

	/**
	 * \brief Return the continuous vector that states this agents move.
	 * 
	 * @return Continuous vector that states this agents move.
	 */
	public ContinuousVector getMovement()
	{
		return _movement;
	}

	/**
	 * \brief Return the index of the grid on which this agent is placed.
	 * 
	 * @return Integer grid index of where this agent is placed.
	 */
	public int getGridIndex()
	{
		return _agentGridIndex;
	}

	/**
	 * \brief Return the LocatedGroup of agents that are present in the
	 * location where this agent is placed.
	 * 
	 * @return	LocatedGroup containing all agents present in the same grid
	 * space as this agent.
	 */
	public LocatedGroup getGridElement()
	{
		return _agentGrid.getShovingGrid()[_agentGridIndex];
	}
	
	/**
	 * \brief Move this agent to another grid index.
	 * 
	 * @param aGridIndex Grid index in which this agent should now be placed.
	 */
	public void setGridIndex(int aGridIndex)
	{
		_agentGridIndex = aGridIndex;
	}

	/**
	 * \brief Return the domain where this agent is contained.
	 * 
	 * @return The domain where this agent is contained (Domain object).
	 */
	public Domain getDomain()
	{
		return _species.domain;
	}
	
	/* bacillus procedures */
	public void randomizeOrientation()
	{
		if (_agentGrid.is3D)
			_headLocation = new ContinuousVector(0d, ExtraMath.random.nextDouble() , ExtraMath.random.nextDouble()); //random orientation; 
		else 
			_headLocation = new ContinuousVector(0.5d , 0.5d, 0d); //random orientation;

		_tailLocation = new ContinuousVector(0d, 0d, 0d);

	}
	
	protected void updateOrientationVector(double angle)
	{
		_headLocation = new ContinuousVector(this._location.x, this._location.y + this._radius - this._capsular_radius, this._location.z); 
		_tailLocation = new ContinuousVector(this._location.x, this._location.y - this._radius + this._capsular_radius, this._location.z);
	
		ContinuousVector[] vo = new ContinuousVector[2];
		vo[0] = this._location;
		vo[1] = new ContinuousVector(this._location.x+1,this._location.y,this._location.z+1);
		
		ContinuousVector[] v = new ContinuousVector[2]; 
		v[0] = this._location;
		
		angle = Math.toRadians(90) ;//Math.random(); //in radians
		
		//rotate head from center
		v[1] = this._headLocation;
		this._headLocation = RotateVector(angle,v,vo);
		
		//rotate tail from center
		v[1] = this._tailLocation;
		this._tailLocation = RotateVector(angle,v,vo);
	}
	
	//return end point
	public static ContinuousVector RotateVector(double theta, 
			ContinuousVector[] v, ContinuousVector[] orientation)
    {
		//0 = start, 1 = end of vector
		//normalization of localized euclidean vector
		double norm = orientation[1].distance(orientation[0]);

		double vo_mag_x = (orientation[1].x - orientation[0].x) / norm;
        double vo_mag_y = (orientation[1].y - orientation[0].y) / norm;
        double vo_mag_z = (orientation[1].z - orientation[0].z) / norm;

        double mag_x = v[1].x-v[0].x;
        double mag_y = v[1].y-v[0].y;
        double mag_z = v[1].z-v[0].z;
        Quaternion Q1 = new Quaternion((double)0,mag_x,mag_y,mag_z);
        
        Quaternion Q2 = new Quaternion((float)Math.cos(theta / 2),
	            (float)(vo_mag_x * Math.sin(theta / 2)),
	            (float)(vo_mag_y * Math.sin(theta / 2)),
	            (float)(vo_mag_z * Math.sin(theta / 2)));
        Quaternion conjQ2 = Quaternion.Conjugate(Q2);

        Quaternion Q3;

        Q3 = Quaternion.Multiply(Quaternion.Multiply(Q2,Q1),conjQ2);

        ContinuousVector result = new ContinuousVector(v[0].x + Q3.x, v[0].y + Q3.y, v[0].z + Q3.z);
        return result;
    }
	

	/**
	 * Apply the rotation angle stored 
	 */
	public void rotate()
	{
	
		if (rotationAngle > 0 && torque.magnitude > 0 && 
				rotationAngle < 0.6)
		{ 
			ContinuousVector[] center_head = {_location,_headLocation};
			ContinuousVector[] center_tail = {_location,_tailLocation};
			ContinuousVector[] T1 = {
	        		new ContinuousVector(torque.start[0],torque.start[1],torque.start[2]),
	        		new ContinuousVector(torque.start[0]+torque.mag_x,
	        				torque.start[1] + torque.mag_y,
	        				torque.start[2] + torque.mag_z)};
	     
		    ContinuousVector newHead = RotateVector(rotationAngle,center_head,T1);
		    ContinuousVector newTail = RotateVector(rotationAngle,center_tail,T1);
		    
		    if (!Double.isNaN(newHead.x) && !Double.isNaN(newTail.x))
		    {
		    	this._headLocation = newHead; 
		    	this._tailLocation = newTail;
		    }
			 
		}
		
		rotationAngle = 0;
		torque = new EuclideanVector(_location,_location);
	}
	
	public void checkBoundariesTailHead() {
		
		ContinuousVector _newTailLoc = new ContinuousVector();
		_newTailLoc.set(_tailLocation);
		_newTailLoc.add(_movement);
		
		ContinuousVector _newHeadLoc = new ContinuousVector();
		_newHeadLoc.set(_headLocation);
		_newHeadLoc.add(_movement);
		
		AllBC aBoundaryT = getDomain().testCrossedBoundary(_newTailLoc);
		AllBC aBoundaryH = getDomain().testCrossedBoundary(_newHeadLoc);
		
		boolean testHead = (aBoundaryH!=null);
		boolean testTail = (aBoundaryT!=null);
		
		
		if (testHead)
		{
			 EuclideanVector force = 
				 new EuclideanVector(_headLocation,aBoundaryH.getOrthoProj(_headLocation));

			 //we need a vector pointing inside with the size of the capsular radius
			 ContinuousVector forceNormal = new ContinuousVector(0d, 0d, 0d); 
			 forceNormal = aBoundaryH.getShape().getNormalInside();
			 forceNormal.normalizeVector();
			 forceNormal.times(_capsular_radius);
			 force.mag_x += forceNormal.x;
			 force.mag_y += forceNormal.y;
			 force.mag_z += forceNormal.z;
			 
			 _movement.add(force.getContinuousVector());
			 
			
			 double[] _center = {_location.x,_location.y,_location.z};
				EuclideanVector N = new  EuclideanVector(force.end,_center);
				EuclideanVector T = force.CrossProduct(N);;
				this.rotationAngle += CollisionEngine.applyForceToCapsule(
						this._location, new EuclideanVector(_tailLocation,_headLocation),
						_capsular_radius, force, -1, null);
				torque = torque.Plus(T);
		}
		
		if (testTail)
		{
			 EuclideanVector force = 
				 new EuclideanVector(_tailLocation,aBoundaryT.getOrthoProj(_tailLocation));
		
		
			 //we need a vector pointing inside with the size of the capsular radius
			 ContinuousVector forceNormal = new ContinuousVector(0d, 0d, 0d); 
			 forceNormal = aBoundaryT.getShape().getNormalInside();
			 forceNormal.normalizeVector();
			 forceNormal.times(_capsular_radius);
			 force.mag_x += forceNormal.x;
			 force.mag_y += forceNormal.y;
			 force.mag_z += forceNormal.z;
			 
			 _movement.add(force.getContinuousVector());
			 
			 double[] _center = {_location.x,_location.y,_location.z};
				EuclideanVector N = new  EuclideanVector(force.end,_center);
				EuclideanVector T = force.CrossProduct(N);;
				this.rotationAngle += CollisionEngine.applyForceToCapsule(
						this._location, new EuclideanVector(_tailLocation,_headLocation),
						_capsular_radius, force, -1, null);
				torque = torque.Plus(T);
		}

	}
	
	public ContinuousVector getHeadLocation() {
		return _headLocation;
	}
	
	public ContinuousVector getTailLocation() {
		return _tailLocation;
	}
}