package be.svlandeg.diffany.usecase.arabidopsis;

import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This data object stores the fold changes, p-values and FDR rates per array ID, for a certain analysis (uniquely identified by name).
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
	 * TODO
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
	 * TODO
	 * @return
	 */
	public SortedSet<String> getArrayIDs()
	{
		return new TreeSet<String>(FDRs.keySet());
	}
	
	/**
	 * TODO
	 * @param arrayID
	 * @return
	 */
	public double getFoldchange(String arrayID)
	{
		return foldchanges.get(arrayID);
	}
	
	/**
	 * TODO
	 * @param arrayID
	 * @return
	 */
	public double getPvalue(String arrayID)
	{
		return pvalues.get(arrayID);
	}
	
	/**
	 * TODO
	 * @param arrayID
	 * @return
	 */
	public double getFDR(String arrayID)
	{
		return FDRs.get(arrayID);
	}
	

}
