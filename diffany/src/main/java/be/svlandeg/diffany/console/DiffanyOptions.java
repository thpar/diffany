package be.svlandeg.diffany.console;

import java.util.HashSet;
import java.util.Set;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 * This class defines the various options available when running the Diffany algorithms through the console.
 * 
 * @author Sofie Van Landeghem
 */
public class DiffanyOptions
{
	
	private Options options;
	
	protected static String helpShort = "h";
	protected static String versionShort = "v";
	protected static String modeShort = "m";
	protected static String logShort = "l";
	
	protected static String refShort = "ref";
	protected static String conShort = "cond";
	protected static String diffShort = "diff";
	protected static String overlapShort = "overlap";

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
		boolean hasArgument = false;
		Set<Option> allFlags = new HashSet<Option>();
		
		allFlags.add(new Option(logShort, "log", hasArgument, "display the log contents"));
		allFlags.add(new Option(versionShort, "version", hasArgument, "print the version information and exit"));
		allFlags.add(new Option(helpShort, "help", hasArgument, "print this help message"));
		
		return allFlags;
	}
	
	/**
	 * Define the options specifying necessary arguments for the Diffany algorithms
	 */
	private Set<Option> getAllParameters()
	{
		boolean hasArgument = true;
		Set<Option> allParameters = new HashSet<Option>();
		
		allParameters.add(new Option(modeShort, "mode", hasArgument, "the mode of the run: pairwise or 1-against-all"));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("reference");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("the input directory containing the reference network");
		allParameters.add(OptionBuilder.create(refShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("conditions");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired(true);
		OptionBuilder.withValueSeparator('|');
		OptionBuilder.withDescription("the input directory containing the condition-specific network(s)");
		allParameters.add(OptionBuilder.create(conShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("differential");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired(true);
		OptionBuilder.withValueSeparator('|');
		OptionBuilder.withDescription("the output directory containing the generated differential network(s)");
		allParameters.add(OptionBuilder.create(diffShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("overlap");
		OptionBuilder.hasArgs();
		OptionBuilder.isRequired(true);
		OptionBuilder.withValueSeparator('|');
		OptionBuilder.withDescription("the output directory containing the generated overlap network(s)");
		allParameters.add(OptionBuilder.create(overlapShort));

		
		
		
		
		return allParameters;
	}
}



