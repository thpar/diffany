package be.svlandeg.diffany.usecase.arabidopsis.osmotic;

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
	 * Retrieve the root directory of the experimental data
	 * @return the root directory, containing subdirectories for the leaf development data and the osmotic stress data
	 */
	public File getInputDataDir()
	{
		return new File(root, "data-marieke"); 
	}
	
	/**
	 * Retrieve the directory containing the leaf development files
	 * @return the directory containing the leaf development files
	 */
	public File getLeafDevelDir()
	{
		return new File(getInputDataDir(), "leaf-development"); 
	}
	
	/**
	 * Retrieve the directory containing the osmotic stress files
	 * @return the directory containing the osmotic stress files
	 */
	public File getRootOsmoticStressDir()
	{
		return new File(getInputDataDir(), "short-term-osmotic-stress"); 
	}
	
}
