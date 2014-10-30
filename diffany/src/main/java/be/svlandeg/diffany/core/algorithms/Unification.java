package be.svlandeg.diffany.core.algorithms;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import be.svlandeg.diffany.core.networks.Edge;
import be.svlandeg.diffany.core.networks.EdgeDefinition;
import be.svlandeg.diffany.core.project.Logger;
import be.svlandeg.diffany.core.semantics.EdgeOntology;

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
	public void expandEdgeOntology(Set<Edge> edges, EdgeOntology eo)
	{
		Map<String, Boolean> theseCats = new HashMap<String, Boolean>();
		Map<String, Boolean> existingCats = new HashMap<String, Boolean>();
		for (Edge e : edges)
		{
			String sourceType = e.getType();
			String sourceCat = eo.getSourceCategory(sourceType);
			
			boolean symmetrical = e.isSymmetrical();
			if (theseCats.containsKey(sourceCat))
			{
				// only symmetrical when all input edges of this type are
				symmetrical = symmetrical && theseCats.get(sourceCat);
			}
			
			// known type
			if (sourceCat != null)
			{
				boolean previousSymmetry = eo.isSymmetricalSourceCat(sourceCat);
				existingCats.put(sourceCat, previousSymmetry);
				theseCats.put(sourceCat, previousSymmetry);
			}
			// unknown type, will be it own class
			else
			{
				theseCats.put(sourceType, symmetrical);
			}
		}
		
		for (String cat : theseCats.keySet())
		{
			boolean symmetryNow = theseCats.get(cat);
			
			// so far unknown type in the ontology
			if (! existingCats.containsKey(cat))
			{
				eo.addSourceCategory(cat, symmetryNow, true);
				eo.addSourceCategoryMapping(cat, cat, false);
				log.log(" Unknown type " + cat + " added to the edge ontology");
			}
			// type already existed in the ontology
			else
			{
				boolean symmetryOld = existingCats.get(cat);
				
				// only symmetrical when all input edges of this type are AND the previous definition is
				boolean symmetryNew = symmetryNow && symmetryOld;
				eo.addSourceCategory(cat, symmetryNew, true);
				if (symmetryNew != symmetryOld)
				{
					log.log(" Edge type " + cat + " changed from directed to symmetrical");
				}
			}
			
		}
		
	}

	/**
	 * Use the directionality of edge types as recorded in the edge ontology to 'fix' a set of edges.
	 * Specifically, symmetrical edges that should be directed, will be split into 2.
	 * 
	 * @param oldEdges the old set of edges
	 * @param eo the edge ontology defining the semantics of edge types
	 * @return the new, unified set of edges adhering to the edge ontology
	 */
	public Set<Edge> unifyEdgeDirection(Set<Edge> oldEdges, EdgeOntology eo)
	{
		Set<Edge> newEdges = new HashSet<Edge>();
		for (Edge oldE : oldEdges)
		{
			EdgeDefinition newDef = new EdgeDefinition(oldE.getDefinition());

			boolean shouldBeSymmetrical = eo.isSymmetricalSourceType(newDef.getType());
			boolean isSymmetrical = newDef.isSymmetrical();
			
			if (isSymmetrical && !shouldBeSymmetrical)
			{
				// source-target: make directed
				newDef.makeSymmetrical(shouldBeSymmetrical);
				
				Edge newEdge1 = new Edge(oldE.getSource(), oldE.getTarget(), newDef);
				newEdges.add(newEdge1);
				
				// make new directed target-source edge
				Edge newEdge2 = new Edge(oldE.getTarget(), oldE.getSource(), newDef);
				newEdges.add(newEdge2);
			}
			else
			{
				// simply keep the current edge
				
				Edge newEdge = new Edge(oldE.getSource(), oldE.getTarget(), newDef);
				newEdges.add(newEdge);
			}
		}
		return newEdges;
	}
}
