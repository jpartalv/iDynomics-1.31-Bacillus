package simulator.geometry.boundaryConditions;

import simulator.SoluteGrid;
import simulator.geometry.Bulk;

public abstract class ConnectedBoundary extends AllBC
{
	/**
	 * The defined bulk in the simulation to which the liquid phase is connected
	 */
	protected Bulk                    _connectedBulk;
	
	/**
	 * \brief Return the bulk that is connected to this boundary
	 * 
	 * TODO Rob 14/1/2015: Generalise this to all computational domains?
	 * 
	 * @return Bulk object that is connected to this boundary
	 */
	public Bulk getBulk()
	{
		return _connectedBulk;
	}
	
	/**
	 * \brief For a specified solute, returns the level of that solute in the
	 * bulk.
	 * 
	 * @param soluteIndex	Index of the solute in the simulation dictionary.
	 * @return	Value of solute in the connected bulk.
	 */
	public Double getBulkValue(int soluteIndex)
	{
		return _connectedBulk.getValue(soluteIndex);
	}
	
	/**
	 * \brief Updates the levels in the bulk.
	 * 
	 *  Allows reaction or flux-based bulk treatment.
	 * 
	 *  @param soluteGrid	Array of all solute grids.
	 *  @param reacGrid		Array of all reaction grids.
	 *  @param timeStep	The internal timestep currently being applied in this 
	 *  simulation.
	 */
	public void updateBulk(SoluteGrid[] soluteGrid,
									SoluteGrid[] reacGrid, Double timeStep)
	{
		_connectedBulk.updateBulk(soluteGrid, reacGrid, timeStep);
	}
}