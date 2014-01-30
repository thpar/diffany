package be.svlandeg.diffany.console;

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
	 * Define the options available for providing meta data (all simple flags)
	 */
	private void defineOptions()
	{
		boolean hasArgument = false;
		
		options = new Options();
		
		options.addOption(new Option(versionShort, "version", hasArgument, "print the version information and exit"));
		options.addOption(new Option(helpShort, "help", hasArgument, "print this help message"));
	}
	
}
