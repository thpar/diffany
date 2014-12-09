package be.svlandeg.diffany.study.osmotic;

import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This data object stores the fold changes, p-values and FDR rates per array ID, for a certain analysis (uniquely identified by name).
 * The data can also be indexed by locus IDs, in which case the boolean value 'arrayIDs' is set to false.
 * 
 * @author Sofie Van Landeghem
 */
public class OverexpressionData
{
	private String name;
	private boolean rawArrayIDs;
	
	private Map<String, Double> foldchanges;
	private Map<String, Double> pvalues;
	private Map<String, Double> FDRs;
	
	/**
	 * Create a new object to contain overexpression values
	 * @param name the (unique) name of this dataset - its uniqueness should be enforced in the use-case
	 * @param rawArrayIDs if true, this data keeps as key raw array IDs (otherwise, they will be locus tags)
	 */
	public OverexpressionData(String name, boolean rawArrayIDs)
	{
		this.name = name;
		this.rawArrayIDs = rawArrayIDs;
		foldchanges = new HashMap<String, Double>();
		pvalues = new HashMap<String, Double>();
		FDRs = new HashMap<String, Double>();
	}
	
	
	/**
	 * Add calculated scores for a certain array ID. This method assumes (and checks) that there were no previous values recorded for this array ID.
	 * If there were, an error message is thrown.
	 * 
	 * @param arrayID the unique array ID (should not have been previously defined)
	 * @param foldchange the fold change
	 * @param pvalue the p-value of the change
	 * @param FDR the false discovery rate of the fold change
	 * 
	 * @throws IllegalArgumentException when values for the arrayID were recorded previously
	 */
	public void addResult(String arrayID, double foldchange, double pvalue, double FDR)
	{
		if (foldchanges.containsKey(arrayID) || pvalues.containsKey(arrayID) || FDRs.containsKey(arrayID))
		{
			String errormsg = "Encountered " + arrayID + " twice in " + name + " ?! ";
			throw new IllegalArgumentException(errormsg);
		}
		
		foldchanges.put(arrayID, foldchange);
		pvalues.put(arrayID, pvalue);
		FDRs.put(arrayID, FDR);
	}
	
	/**
	 * Return whether or not this dataset is indexed by array IDs (if not, they are locus IDs)
	 * @return whether or not this dataset is indexed by array IDs
	 */
	public boolean indexedByRawArrayIDs()
	{
		return rawArrayIDs;
	}
	
	/**
	 * Return the name of this data
	 * @return the name
	 */
	public String getName()
	{
		return name;
	}
	
	/**
	 * Retrieve all IDs indexed in this dataset
	 * @return all indexed IDs
	 */
	public SortedSet<String> getArrayIDs()
	{
		return new TreeSet<String>(FDRs.keySet());
	}
	
	/**
	 * Retrieve the fold change of a certain ID
	 * @param arrayID the query ID
	 * @return the fold change (FC) of this ID
	 */
	public double getFoldchange(String arrayID)
	{
		return foldchanges.get(arrayID);
	}
	
	/**
	 * Retrieve the p-value of a certain ID
	 * @param arrayID the query ID
	 * @return the p-value of this ID
	 */
	public double getPvalue(String arrayID)
	{
		return pvalues.get(arrayID);
	}
	
	/**
	 * Retrieve the false discovery rate (FDR) of a certain ID
	 * @param arrayID the query ID
	 * @return the false discovery rate (sort of normalized p-value) of this ID
	 */
	public double getFDR(String arrayID)
	{
		return FDRs.get(arrayID);
	}
	
	/**
	 * Retrieve all the significant gene IDs in this overexpression dataset, by using the threshold as a minimal cutoff of the FDR values.
	 *
	 * @param threshold the FDR cutoff
	 * @return all node IDs above the threshold, mapped to their corresponding fold change
	 */
	public Map<String, Double> getSignificantGenes(double threshold)
	{
		Map<String, Double> nodes = new HashMap<String, Double>();

		for (String id : getArrayIDs())
		{
			double FDR = getFDR(id);
			if (FDR <= threshold)
			{
				double FC = getFoldchange(id);
				nodes.put(id.toLowerCase(),  FC);
			}
		}
		return nodes;
	}
	

}
