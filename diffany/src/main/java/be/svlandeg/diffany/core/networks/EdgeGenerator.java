package be.svlandeg.diffany.core.networks;

/**
 * This class generates special edges to be used in the input/output networks, such as default edges,
 * virtual edges or void edges.
 * 
 * @author Sofie Van Landeghem
 */
public class EdgeGenerator
{
	
	protected static String VOID_DIRECT_TYPE = "*nodirectedtype*";
	protected static String VOID_SYMM_TYPE = "*nosymmetricaltype*";
	
	
	/**
	 * Obtain a default edge. It will be initalized to default values from {@link EdgeDefinition}: "unspecified_connection", weight 1, symmetrical, not negated.
	 * 
	 * @return a default edge (weight == 1, symmetrical == true)
	 */
	public EdgeDefinition getDefaultEdge()
	{
		return new EdgeDefinition(EdgeDefinition.DEFAULT_TYPE, EdgeDefinition.DEFAULT_SYMM, EdgeDefinition.DEFAULT_WEIGHT, EdgeDefinition.DEFAULT_NEG);
	}
	
	
	/**
	 * Obtain a void edge, for the purpose of being able to compare it to existing edges.
	 * 
	 * @param symmetrical whether or not the edge should be symmetrical
	 * @return a void edge (weight == 0)
	 */
	public EdgeDefinition getVoidEdge(boolean symmetrical)
	{
		if (symmetrical)
		{
			return new EdgeDefinition(VOID_SYMM_TYPE , true, 0, EdgeDefinition.DEFAULT_NEG);
		}
		return new EdgeDefinition(VOID_DIRECT_TYPE , false, 0, EdgeDefinition.DEFAULT_NEG);
	}
	
	
	

}
