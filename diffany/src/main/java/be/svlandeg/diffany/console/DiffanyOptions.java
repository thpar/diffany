package be.svlandeg.diffany.console;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 * This class defines the various options available when running the Diffany algorithms through the console.
 * 
 * @author Sofie Van Landeghem
 */
public class DiffanyOptions
{
	
	private Options options;
	protected static String logShort = "l";
	protected static String refShort = "ref";

	/**
	 * Constructor initializes the options available in Diffany
	 */
	public DiffanyOptions()
	{
		defineOptions();
	}

	/**
	 * Retrieve the options object for the Diffany project
	 * @return the Diffany options object
	 */
	public Options getDiffanyOptions()
	{
		return options;

	}
	
	/**
	 * Define the options available in the Diffany project
	 */
	private void defineOptions()
	{
		options = new Options();
		for (Option o : getAllFlags())
		{
			options.addOption(o);
		}
		for (Option o : getAllParameters())
		{
			options.addOption(o);
		}
	}
	
	/**
	 * Define the flag options available in the Diffany project
	 */
	private Set<Option> getAllFlags()
	{
		boolean argument = false;
		Set<Option> allFlags = new HashSet<Option>();
		
		allFlags.add(new Option(logShort, "log", argument, "display the log contents"));
		
		return allFlags;
	}
	
	/**
	 * Define the options specifying necessary arguments for the Diffany algorithms
	 */
	private Set<Option> getAllParameters()
	{
		boolean argument = true;
		Set<Option> allParameters = new HashSet<Option>();
		
		allParameters.add(new Option(refShort, "reference", argument, "the reference network "));
		
		return allParameters;
	}
}



