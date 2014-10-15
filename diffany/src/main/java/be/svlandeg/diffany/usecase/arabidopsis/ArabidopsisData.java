package be.svlandeg.diffany.usecase.arabidopsis;

import java.net.URI;
import java.net.URISyntaxException;

/**
 * This class defines and analyses the data retrieved for Arabidopsis thaliana.
 * 
 * @author Sofie Van Landeghem
 */
public class ArabidopsisData
{
	
	private static String cornetPPIFile = "validated_cornet_all_ppi_table_17012012.tab";
	private static String cornetRegFile = "reg_net_20100205.tab";
	private static String kinaseInteractionFile = "kinase-targets_20131210.csv";
	private static String kinaseFunctionFile = "kinase_activity_go_15102014.tab";
	private static String phosphatFile = "phosphat_20130429.csv";
	

	/**
	 * Retrieve the URI of the CORNET PPI data
	 * @return the URI of the PPI data, or null if the resource could not be located
	 */
	public URI getCornetPPI()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + cornetPPIFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + cornetPPIFile);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the CORNET regulatory data
	 * @return the URI of the regulatory data, or null if the resource could not be located
	 */
	public URI getCornetReg()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + cornetRegFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + cornetRegFile);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the Phosphat data
	 * @return the URI of the Phosphat data, or null if the resource could not be located
	 */
	public URI getPhosphat()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + phosphatFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + phosphatFile);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the kinase activity data
	 * @return the URI of the kinase activity data, or null if the resource could not be located
	 */
	public URI getKinases()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + kinaseFunctionFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + kinaseFunctionFile);
        }
		return null;
	}
	
	/**
	 * Retrieve the URI of the kinase interaction data
	 * @return the URI of the kinase interaction data, or null if the resource could not be located
	 */
	public URI getKinaseInteractions()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + kinaseInteractionFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + kinaseInteractionFile);
        }
		return null;
	}
}
