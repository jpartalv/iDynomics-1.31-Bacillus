package simulator.agent.zoo;

import java.util.ArrayList;
import java.util.LinkedList;

import simulator.Simulator;
import simulator.agent.ActiveParam;
import simulator.agent.zoo.PlasmidBac;
import simulator.agent.Species;
import utils.LogFile;
import utils.XMLParser;

/**
 * \brief TODO
 * 
 * @author Robert Clegg (r.j.clegg@bham.ac.uk)
 */
public class PlasmidParam extends ActiveParam
{
	/**
	 * Info used in output writing.
	 */
	public String name;
	
	/**
	 * Default copy number of this plasmid in a host.
	 */
	public Integer copyNumDefault = 1;
	
	/**
	 * Length of the pilus associated with this plasmid for conjugation 
	 * (in um). Note that the pilus should reach between cell surfaces, not
	 * between cell centres.
	 */
	public Double pilusLength = 5.0;
	
	/**
	 * After donating a plasmid of this species, a host needs to recover
	 * before it may exchange again. This parameter gives the time needed to
	 * recover (in hours).
	 * 
	 * TODO Note that this was previously called transferLag - protocol files
	 * will need to be modified.
	 */
	public Double donationLag = 5.0;
	
	/**
	 * After receiving a plasmid of this species, a host needs to recover
	 * before it may receive again. This parameter gives the time needed to
	 * recover (in hours).
	 */
	public Double receptionLag = 1.0;
	
	/**
	 * Probability that this plasmid will be lost by the host during division.
	 */
	public Double lossProbability = 0.0;
	
	/**
	 * Probability that a possible transfer will be successful.
	 */
	public Double transferProficiency = 1.0;
	
	/**
	 * The maximum scan speed (in cells/hour) of Plasmids of this species.
	 */
	public Double scanSpeed = 10.0;
	
	/**
	 * Plasmid markers (used in host compatibility).
	 */
	public ArrayList<String> hostCompatibilityMarkers =
													new ArrayList<String>();
	
	/**
	 * Plasmid Markers for plasmid compatibility (based on incompatibility
	 * groups).
	 */
	public ArrayList<String> plasmidCompatibilityMarkers =
													new ArrayList<String>();
	
	/**
	 * List of integer indices corresponding to the reactions that Plasmids
	 * of this species encode for.
	 */
	public ArrayList<Integer> reactionsEncoded = new ArrayList<Integer>();
	
	
	/*************************************************************************
	 * CONSTRUCTORS
	 ************************************************************************/
	
	public PlasmidParam()
	{
		super();
	}
	
	/**
	 * Called during creation of the species
	 */
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim, aSpeciesRoot, speciesDefaults);
		
		Integer tempInt;
		Double tempDbl;
		LinkedList<XMLParser> tempXML;
		
		name = aSpeciesRoot.getName();
		
		tempInt = getSpeciesParameterInteger("copyNumDefault",
											aSpeciesRoot, speciesDefaults);
		copyNumDefault = ( tempInt == null ) ? copyNumDefault : tempInt;
		
		tempDbl = getSpeciesParameterLength("pilusLength", 
											aSpeciesRoot, speciesDefaults);
		pilusLength = Double.isFinite(tempDbl) ? tempDbl : pilusLength;
		
		tempDbl = getSpeciesParameterTime("donationLag",
											aSpeciesRoot, speciesDefaults);
		donationLag = Double.isFinite(tempDbl) ? tempDbl : donationLag;
		
		tempDbl = getSpeciesParameterTime("receptionLag",
											aSpeciesRoot, speciesDefaults);
		receptionLag = Double.isFinite(tempDbl) ? tempDbl : receptionLag;
		
		tempDbl = getSpeciesParameterDouble("lossProbability",
											aSpeciesRoot, speciesDefaults);
		lossProbability = Double.isFinite(tempDbl)? tempDbl : lossProbability;

		tempDbl = getSpeciesParameterDouble("transferProficiency",
											aSpeciesRoot, speciesDefaults);
		transferProficiency = Double.isFinite(tempDbl) ?
												tempDbl : transferProficiency;
		tempDbl = getSpeciesParameterTime("scanSpeed",
											aSpeciesRoot, speciesDefaults);
		scanSpeed = Double.isFinite(tempDbl) ? tempDbl : scanSpeed;
		/*
		 * 
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("reaction");
		for ( XMLParser parser : tempXML )
			if ( parser.getAttribute("status").equals("active") )
				reactionsEncoded.add(aSim.getReactionIndex(parser.getName()));
		/*
		 * If no host markers are given, assume all species may be hosts.
		 * Note: Old EpiBac protocol files should be unaffected by this.
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("compatibleHost");
		if ( tempXML == null || tempXML.isEmpty() )
		{
			//System.out.println("Entering all host compatability markers");
			for ( Species spec : aSim.speciesList )
				if ( spec.getProgenitor() instanceof PlasmidBac )
					hostCompatibilityMarkers.add( spec.speciesName );
		}
		else
		{
			//System.out.println("Trying to read in host compatability markers");
			for ( XMLParser parser : tempXML )
			{
				//System.out.println("\tAdding "+parser.getName());
				hostCompatibilityMarkers.add( parser.getName() );
			}
		}
		/*
		 * Try to add this plasmid name to the list of potential plasmids in 
		 * compatible hosts, if they have been initialised yet.
		 */
		for ( String hostName : hostCompatibilityMarkers )
			try
			{
				((PlasmidBacParam)aSim.getSpecies(hostName).getSpeciesParam())
											.potentialPlasmids.add(this.name);
			}
			catch ( Exception e) {}
		/*
		 * If no plasmid markers are given, assume all plasmids are
		 * incompatible.
		 * 
		 * TODO check that compatibility goes both ways?
		 */
		tempXML = aSpeciesRoot.getChildrenParsers("compatiblePlasmids");
		if ( ! (tempXML == null || tempXML.isEmpty()) )
			for ( XMLParser parser : tempXML )
				plasmidCompatibilityMarkers.add( parser.getName() );
	}
	
	/**
	 * \brief Check if a potential host is compatible with a plasmid of this
	 * species.
	 * 
	 * @param targetRecipient PlasmidBac to check compatibility with.
	 * @return boolean: true if target is compatible, false if not.
	 */
	public boolean isCompatible(PlasmidBac targetRecipient)
	{
		if ( ! hostCompatibilityMarkers.contains(targetRecipient.getName()) )
			return false;
		for ( Plasmid p : targetRecipient.getPlasmidsHosted() )
			if ( ! plasmidCompatibilityMarkers.contains(p.getName()) )
				return false;
		return true;
	}
}
