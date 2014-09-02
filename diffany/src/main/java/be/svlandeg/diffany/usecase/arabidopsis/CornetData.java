package be.svlandeg.diffany.usecase.arabidopsis;

import java.net.URI;
import java.net.URISyntaxException;

/**
 * This class defines and analyses the data retrieved from CORNET.
 * 
 * @author Sofie Van Landeghem
 */
public class CornetData
{
	
	private static String cornetPPIDataFile = "validated_cornet_all_ppi_table_17012012.tab";
	private static String cornetRegDataFile = "reg_net_20100205.tab";

	/**
	 * Retrieve the URI of the CORNET PPI data
	 * @return the URI of the PPI data, or null if the resource could not be located
	 */
	public URI getCornetPPI()
	{
		try
        {
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + cornetPPIDataFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + cornetPPIDataFile);
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
	        return Thread.currentThread().getContextClassLoader().getResource("data/" + cornetRegDataFile).toURI();
        }
        catch (URISyntaxException e)
        {
	        System.out.println(" !  Couldn't read " + cornetRegDataFile);
        }
		return null;
	}
}
