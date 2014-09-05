package be.svlandeg.diffany.usecase;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;

import be.svlandeg.diffany.core.networks.Node;
import be.svlandeg.diffany.usecase.arabidopsis.GenePrinter;
import be.svlandeg.diffany.usecase.arabidopsis.OverexpressionData;


/**
 * This class is used to analyse expression datasets in the Diffany format.
 * 
 * @author Sofie Van Landeghem
 */
public class ExpressionDataAnalysis
{

	
	/**
	 * Constructor: currently empty
	 */
	public ExpressionDataAnalysis()
	{}
	
	/**
	 * Retrieve all the significant gene IDs in an overexpression dataset, by using the threshold as a minimal cutoff of the FDR values.
	 * @param data the input datasets
	 * @param threshold the FDR cutoff
	 * @return all node IDs above the threshold, mapped to their corresponding fold change
	 */
	public Map<String, Double> getSignificantGenes(OverexpressionData data, double threshold)
	{
		Map<String, Double> nodes = new HashMap<String, Double>();

		SortedSet<String> ids = data.getArrayIDs();
		for (String id : ids)
		{
			double FDR = data.getFDR(id);
			if (FDR <= threshold)
			{
				double FC = data.getFoldchange(id);
				nodes.put(id.toLowerCase(),  FC);
			}
		}
		return nodes;
	}

}
