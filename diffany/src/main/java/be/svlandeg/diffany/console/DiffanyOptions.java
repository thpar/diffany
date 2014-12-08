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
	protected static String condShort = "cond";
	protected static String outputShort = "output";
	
	protected static String diffNameShort = "diffName";
	protected static String consNameShort = "consName";
	
	protected static String diffID = "diffID";
	protected static String consID = "consID";
	
	protected static String cutoffShort = "conf";
	protected static String operatorShort = "oper";
	protected static String modeShort = "mode";
	
	protected static boolean defaultMinOperator = true;
	protected static boolean defaultModePairwise = false;
	

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
		
		allFlags.add(new Option(logShort, "log", hasArgument, "display a progress/log file during the run"));
		
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
		OptionBuilder.withLongOpt("referenceDir");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("the input directory containing the reference network");
		allParameters.add(OptionBuilder.create(refShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("conditionsDir");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired(true);
		//OptionBuilder.withValueSeparator('|');
		OptionBuilder.withDescription("the input directory containing the condition-specific network");
		allParameters.add(OptionBuilder.create(condShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("outputDir");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired(true);
		OptionBuilder.withDescription("the output directory which will contain the generated differential/consensus networks");
		allParameters.add(OptionBuilder.create(outputShort));

		OptionBuilder.withLongOpt("differentialName");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the name of the generated differential network");
		allParameters.add(OptionBuilder.create(diffNameShort));
		
		OptionBuilder.withLongOpt("consensusName");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the name of the generated consensus network");
		allParameters.add(OptionBuilder.create(consNameShort));
		
		OptionBuilder.withLongOpt("differentialID");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the ID of the differential network - it will not be generated if this is a negative value");
		allParameters.add(OptionBuilder.create(diffID));
		
		OptionBuilder.withLongOpt("consensusID");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the ID of the consensus network - it will not be generated if this is a negative value");
		allParameters.add(OptionBuilder.create(consID));
		
		OptionBuilder.withLongOpt("confidenceMin");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the minimum confidence threshold for differential and consensus edges, as an integer or double (default=0.0)");
		allParameters.add(OptionBuilder.create(cutoffShort));
		
		String defaultMinOperatorString = "min";
		if (! defaultMinOperator)
		{
			defaultMinOperatorString = "max";
		}
		
		OptionBuilder.withLongOpt("operator");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the operator used to create consensus edges: min or max (default=" + defaultMinOperatorString + ")");
		allParameters.add(OptionBuilder.create(operatorShort));
		
		String defaultModeString = "pairwise";
		if (! defaultModePairwise)
		{
			defaultMinOperatorString = "all";
		}
		
		OptionBuilder.withLongOpt("comparisonMode");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the mode of comparison: pairwise or all (default=" + defaultModeString + ")");
		allParameters.add(OptionBuilder.create(modeShort));
		
		return allParameters;
	}
}



