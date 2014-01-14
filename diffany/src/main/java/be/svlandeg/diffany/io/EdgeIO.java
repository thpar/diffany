package be.svlandeg.diffany.io;

import java.text.DecimalFormat;

import be.svlandeg.diffany.concepts.Edge;
import be.svlandeg.diffany.concepts.EdgeDefinition;

/**
 * This class allows reading or writing an {@link Edge} from File.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeIO
{
	

	/**
	 * Get a string representation of an edge.
	 * More specifically, print it as: source.name - target.name - edge.type - symmetrical - weight - negated.
	 * 
	 * @return a string representation of this edge, ready for printing
	 */
	public static String writeToTab(Edge e)
	{
		String defResult = writeDefinitionToTab(e);
		String result = e.getSource().getName() + '\t' + e.getTarget().getName() + '\t' + defResult;
		return result;
	}
	
	/**
	 * Get a string representation of an edge definition.
	 * More specifically, print it as: edge.type - symmetrical - weight - negated.
	 * 
	 * @return a string representation of this edge, ready for printing
	 */
	public static String writeDefinitionToTab(EdgeDefinition def)
	{
		DecimalFormat df = new DecimalFormat("#.##");
		String symm = "symmetrical";
		if (! def.isSymmetrical())
		{
			symm = "directed";
		}
		String neg = "negated";
		if (! def.isNegated())
		{
			neg = "not negated";
		}
		double weight = def.getWeight();
		String result = def.getType() + '\t' + symm + '\t' + df.format(weight) + '\t' + neg;
		
		return result;
	}
	
	/**
	 * Read an Edge from a tab-delimited String.
	 * 
	 * @param s the tab-delimited string containing all parameters of the Edge
	 * @return the corresponding Edge object
	 */
	public Edge readFromTab(String s)
	{
		//TODO v1.1: implement!
		throw new UnsupportedOperationException("readFromTab not yet implemented");
	}
	
	

}
