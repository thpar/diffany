package be.svlandeg.diffany.usecase.osmotic;

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
	 * Retrieve the directory containing the leaf development files
	 * @return
	 */
	public File getInputDataDir()
	{
		return new File(root, "data-marieke"); 
	}
	
	/**
	 * Retrieve the directory containing the leaf development files
	 * @return
	 */
	public File getLeafDevelDir()
	{
		return new File(getInputDataDir(), "leaf-development"); 
	}
	
	/**
	 * Retrieve the directory containing the osmotic stress files
	 * @return
	 */
	public File getOsmoticStressDir()
	{
		return new File(getInputDataDir(), "short-term-osmotic-stress"); 
	}
}
