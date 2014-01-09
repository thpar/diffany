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
		this.log = log;
	}
	
	/**
	 * Expand an existing edge ontology to cover all interaction types from the input data
	 * @param edges the set of edges in the input network
	 * @param eo the edge ontology
	 */
	public void expandOntology(Set<Edge> edges, EdgeOntology eo)
	{
		Map<String, Boolean> newTypes = new HashMap<String, Boolean>();
		for (Edge e : edges)
		{
			String sourceType = e.getType();
			String sourceClass = eo.getSourceCategory(sourceType);
			if (sourceClass == null)
			{
				boolean symmetrical = e.isSymmetrical();
				if (newTypes.containsKey(sourceType))
				{
					// only symmetrical when all input edges of this type are
					symmetrical = symmetrical && newTypes.get(sourceType);
				}
				newTypes.put(sourceType, symmetrical);
			}
		}
		
		for (String type : newTypes.keySet())
		{
			eo.addSourceCategory(type, newTypes.get(type));
			eo.addSourceCategoryMapping(type, type, false);
			log.log(" Unknown type " + type + " added to the edge ontology.");
		}
	}
}
