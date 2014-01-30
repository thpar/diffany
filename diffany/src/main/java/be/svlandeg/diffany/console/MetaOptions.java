package be.svlandeg.diffany.console;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

/**
 * This class defines the various 'meta' options available when running the Diffany jar through the console,
 * such as the help message and the version number.
 * 
 * @author Sofie Van Landeghem
 */
public class MetaOptions
{

	private Options options;
	
	protected static String helpShort = "h";
	protected static String versionShort = "v";

	/**
	 * Constructor initializes the meta options for this program
	 */
	public MetaOptions()
	{
		defineOptions();
	}
	
	/**
	 * Retrieve the options object for the meta data
	 * @return the options object
	 */
	public Options getMetaOptions()
	{
		return options;

	}

	
	/**
	 * Define the options available for providing meta data
	 */
	private void defineOptions()
	{
		options = new Options();
		for (Option o : getAllFlags())
		{
			options.addOption(o);
		}
	}
	
	/**
	 * Define the flag options available in the Diffany project
	 */
	private Set<Option> getAllFlags()
	{
		boolean hasArgument = false;
		Set<Option> allFlags = new HashSet<Option>();
		
		allFlags.add(new Option(versionShort, "version", hasArgument, "print the version information and exit"));
		allFlags.add(new Option(helpShort, "help", hasArgument, "print this help message"));
		
		return allFlags;
	}
	
}
