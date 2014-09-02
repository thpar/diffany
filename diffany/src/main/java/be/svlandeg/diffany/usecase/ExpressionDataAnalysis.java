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
	
	private GenePrinter gp;
	
	/**
	 * Constructor: defines the gene printer that can deal with gene synonymy etc.
	 * @param gp the gene printer object
	 */
	public ExpressionDataAnalysis(GenePrinter gp)
	{
		this.gp = gp;
	}
	
	/**
	 * Retrieve all the significant genes in an overexpression dataset, by using the threshold as a minimal cutoff of the FDR values.
	 * @param data the input datasets
	 * @param threshold the FDR cutoff
	 * @return all nodes above the threshold, mapped to their corresponding fold change
	 */
	public Map<Node, Double> getSignificantGenes(OverexpressionData data, double threshold)
	{
		boolean arrayID = data.indexedByRawArrayIDs();

		Map<Node, Double> nodes = new HashMap<Node, Double>();

		SortedSet<String> ids = data.getArrayIDs();
		for (String id : ids)
		{
			double FDR = data.getFDR(id);
			if (FDR <= threshold)
			{
				String symbol = gp.getSymbolByLocusID(id);
				if (arrayID)
				{
					symbol = Arrays.toString(gp.getSymbolByArrayID(id).toArray());
				}
				if (symbol == null)
				{
					symbol = id;
				}
				double FC = data.getFoldchange(id);
				nodes.put(new Node(id.toLowerCase(), symbol, false), FC);
			}
		}
		return nodes;
	}

}
