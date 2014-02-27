package be.svlandeg.diffany.algorithms;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.Logger;
import be.svlandeg.diffany.semantics.EdgeOntology;

/**
 * This class provides generic methods useful for network unification before or
 * after applying differential algorithms.
 * 
 * @author Sofie Van Landeghem
 */
public class Unification
{
	
	private Logger log;
	
	/**
	 * Create a new unification object, which can log important messages.
	 * @param log the logger object
	 */
	public Unification(Logger log)
	{
		if (log == null)
		{
			String errormsg = "The logger should not be null!";
			throw new IllegalArgumentException(errormsg);
		}
		this.log = log;
	}
	
	/**
	 * Expand an existing edge ontology to cover all interaction types from the input data.
	 * Adjust the directionality when need be (e.g. when seeing directed edges for a type that was thought to be symmetrical).
	 * Important events are logged.
	 * 
	 * @param edges the set of edges in the input network
	 * @param eo the edge ontology
	 */
	public void expandOntology(Set<Edge> edges, EdgeOntology eo)
	{
		Map<String, Boolean> theseTypes = new HashMap<String, Boolean>();
		Map<String, Boolean> existingTypes = new HashMap<String, Boolean>();
		for (Edge e : edges)
		{
			String sourceType = e.getType();
			String sourceClass = eo.getSourceCategory(sourceType);
			
			boolean symmetrical = e.isSymmetrical();
			if (theseTypes.containsKey(sourceType))
			{
				// only symmetrical when all input edges of this type are
				symmetrical = symmetrical && theseTypes.get(sourceType);
			}
			theseTypes.put(sourceType, symmetrical);
			
			// known type
			if (sourceClass != null)
			{
				boolean previousSymmetry = eo.isSymmetricalSourceType(sourceType);
				existingTypes.put(sourceType, previousSymmetry);
			}
		}
		
		for (String type : theseTypes.keySet())
		{
			boolean symmetryNow = theseTypes.get(type);
			
			// so far unknown type in the ontology
			if (! existingTypes.containsKey(type))
			{
				eo.addSourceCategory(type, symmetryNow, true);
				eo.addSourceCategoryMapping(type, type, false);
				log.log(" Unknown type " + type + " added to the edge ontology");
			}
			// type already existed in the ontology
			else
			{
				boolean symmetryOld = existingTypes.get(type);
				
				// only symmetrical when all input edges of this type are ánd the previous definition is
				boolean symmetryNew = symmetryNow && symmetryOld;
				eo.addSourceCategory(type, symmetryNew, true);
				if (symmetryNew != symmetryOld)
				{
					log.log(" Edge type " + type + " changed from directed to symmetrical");
				}
			}
			
		}
		
	}
}
