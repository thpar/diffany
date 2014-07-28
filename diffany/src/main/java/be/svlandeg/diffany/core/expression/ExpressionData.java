package be.svlandeg.diffany.core.expression;

import java.util.List;

/**
 * This class models an expression dataset
 * 
 * @author Sofie Van Landeghem
 */
public class ExpressionData
{
	
	private String collectionName;
	private List<String> rows_genes;
	private List<String> columns_samples;
	private double[][] expvalues;
	private boolean normalized;
	
	/**
	 * Create a new expression dataset
	 * @param collectionName the name of the experiment
	 * @param rows_genes the names of the rows (genes)
	 * @param columns_samples the names of the columns (samples)
	 * @param expvalues the actual expression values
	 * @param normalized whether or not the data is already normalized
	 */
	public ExpressionData(String collectionName, List<String> rows_genes, List<String> columns_samples, double[][] expvalues, boolean normalized)
	{
		if (rows_genes.size() != expvalues.length)
		{
			String errormsg = "The number of rows does not match the number of gene names";
			throw new IllegalArgumentException(errormsg);
		}
		if (columns_samples.size() != expvalues[0].length)
		{
			String errormsg = "The number of columns does not match the number of samples";
			throw new IllegalArgumentException(errormsg);
		}
		this.collectionName = collectionName;
		this.rows_genes = rows_genes;
		this.columns_samples = columns_samples;
		this.expvalues = expvalues;
		this.normalized = normalized;
	}

	/**
	 * Return the collection name of this experiment
	 * @return the collection name of this experiment
	 */
	public String getCollectionName()
	{
		return collectionName;
	}

	/**
	 * Return the names/identifiers of the genes (rows)
	 * @return the gene names
	 */
	public List<String> getGenes()
	{
		return rows_genes;
	}

	/**
	 * Return the names of the samples (columns)
	 * @return the sample names
	 */
	public List<String> getSamples()
	{
		return columns_samples;
	}

	/**
	 * Return the actual expression values
	 * @return the expression values
	 */
	public double[][] getExpvalues()
	{
		return expvalues;
	}

	/**
	 * Check whether the expression values can be considered to be normalized already or not
	 * @return whether or not the data is already normalized
	 */
	public boolean isNormalized()
	{
		return normalized;
	}
	
	public String toString()
	{
		String resultString = "Expression data for experiment " + getCollectionName();
		resultString += " - " + getGenes().size() + " genes";
		resultString += " - " + getSamples().size() + " samples";
		return resultString;
	}

}
