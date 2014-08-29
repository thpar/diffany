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
	
	//protected static String modeShort = "m";
	protected static String logShort = "l";
	
	protected static String refShort = "ref";
	protected static String conShort = "cond";
	protected static String diffShort = "diff";
	protected static String consensusShort = "consensus";

	protected static String diffnameShort = "name";
	protected static String cutoffShort = "conf";
	
	protected static String diffID = "diffID";
	protected static String consensusID = "consensusID";

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
		
		allFlags.add(new Option(logShort, "log", hasArgument, "display a log file after running the algorithm"));
		
		return allFlags;
	}
	
	/**
	 * Define the options specifying necessary arguments for the Diffany algorithms
	 */
	private Set<Option> getAllParameters()
	{
		Set<Option> allParameters = new HashSet<Option>();
		
		//allParameters.add(new Option(modeShort, "mode", hasArgument, "the mode of the run: pairwise or 1-against-all"));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("referenceDirectory");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("the input directory containing the reference network");
		allParameters.add(OptionBuilder.create(refShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("conditionsDirectory");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired(true);
		//OptionBuilder.withValueSeparator('|');
		OptionBuilder.withDescription("the input directory containing the condition-specific network");
		allParameters.add(OptionBuilder.create(conShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("differentialDirectory");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired(true);
		OptionBuilder.withDescription("the output directory which will contain the generated differential network");
		allParameters.add(OptionBuilder.create(diffShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("consensusDirectory");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired(true);
		OptionBuilder.withDescription("the output directory which will contain the generated consensus network");
		allParameters.add(OptionBuilder.create(consensusShort));

		OptionBuilder.withLongOpt("networkOutputName");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the name of the generated differential network");
		allParameters.add(OptionBuilder.create(diffnameShort));
		
		OptionBuilder.withLongOpt("diffID");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the ID of the generated differential network");
		allParameters.add(OptionBuilder.create(diffID));
		
		OptionBuilder.withLongOpt("consensusID");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the ID of the generated consensus network");
		allParameters.add(OptionBuilder.create(consensusID));
		
		OptionBuilder.withLongOpt("confidenceMin");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the minimum confidence threshold for differential and consensus edges");
		allParameters.add(OptionBuilder.create(cutoffShort));
		
		return allParameters;
	}
}



