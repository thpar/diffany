package be.svlandeg.diffany.usecase.arabidopsis.tf;

import java.io.File;


/**
 * Class that handles IO for this specific use-case, assuming a fixed directory structure
 * when the root location is known.
 * 
 * @author Sofie Van Landeghem
 */
public class DataIO
{
	
	private String root;
	
	/**
	 * Create a new DataIO object, specifying the root directory which contains the 
	 * leaf development and osmotic stress raw data files
	 * 
	 * @param root the root directory of the input data
	 */
	public DataIO(String root)
	{
		this.root = root;
	}
	
	/**
	 * Retrieve the root directory of the TF-target data
	 * @return the root directory, containing the input TF files
	 */
	public File getTFInputDataDir()
	{
		return new File(root, "data-jan"); 
	}
	
	/**
	 * Retrieve the root directory of the expression data
	 * @return the root directory, containing the input expression files
	 */
	public File getExpInputDataDir()
	{
		return new File(root, "CORNET2.0"); 
	}
	
	/**
	 * Retrieve a file with +/- 1100 experimentally validated TF-target interactions
	 * @return the root directory, containing the input files
	 */
	public File getTFs()
	{
		return new File(getTFInputDataDir(), "AGRIS_CELLWALL.txt"); 
	}

}
