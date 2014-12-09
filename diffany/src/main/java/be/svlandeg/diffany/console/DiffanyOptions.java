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
	
	protected static String logShort = "l";
	
	protected static String cutoffShort = "c";
	protected static String operatorShort = "oper";
	protected static String modeShort = "m";
	
	protected static String headerShort = "h";
	
	protected static String inputShort = "i";
	protected static String outputShort = "o";
	
	protected static String runDiff = "diff";
	protected static String runCons = "cons";
	protected static String nextID = "ID";
	
	protected static boolean defaultRunDiff = true;
	protected static boolean defaultRunCons = true;
	
	protected static boolean defaultMinOperator = true;
	protected static boolean defaultModePairwise = false;
	
	protected static boolean defaultReadHeader = true;
	

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
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("inputDir");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("the input directory containing the reference and condition-specific networks");
		allParameters.add(OptionBuilder.create(inputShort));
		
		OptionBuilder.withArgName("dir");
		OptionBuilder.withLongOpt("outputDir");
		OptionBuilder.hasArgs(1);
		OptionBuilder.isRequired();
		OptionBuilder.withDescription("the output directory which will contain the generated differential/consensus networks");
		allParameters.add(OptionBuilder.create(outputShort));
		
		String defaultRunDiffString = defaultRunDiff? "yes" : "no";
		OptionBuilder.withLongOpt("differential");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("whether or not to calculate differential networks: yes or no (default=" + defaultRunDiffString + ")");
		allParameters.add(OptionBuilder.create(runDiff));
		
		String defaultRunConsString = defaultRunCons? "yes" : "no"; ;
		OptionBuilder.withLongOpt("consensus");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("whether or not to calculate consensus networks: yes or no (default=" + defaultRunConsString + ")");
		allParameters.add(OptionBuilder.create(runCons));
		
		OptionBuilder.withLongOpt("outputID");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the first ID that will be used for the generated networks");
		allParameters.add(OptionBuilder.create(nextID));
		
		OptionBuilder.withLongOpt("confidence");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the minimum confidence threshold for output edges, as an integer or double (default=0.0)");
		allParameters.add(OptionBuilder.create(cutoffShort));
		
		String defaultMinOperatorString = defaultMinOperator? "min" : "max";
		OptionBuilder.withLongOpt("operator");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the operator used to create consensus edges: min or max (default=" + defaultMinOperatorString + ")");
		allParameters.add(OptionBuilder.create(operatorShort));
		
		String defaultModeString = defaultModePairwise? "pairwise" : "all";
		OptionBuilder.withLongOpt("mode");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("the mode of comparison: pairwise or all (default=" + defaultModeString + ")");
		allParameters.add(OptionBuilder.create(modeShort));
		
		String defaultHeaderString = defaultReadHeader? "yes" : "no";
		OptionBuilder.withLongOpt("skipHeader");
		OptionBuilder.hasArgs(1);
		OptionBuilder.withDescription("whether or not to skip the first line (header) in the network .txt files (default=" + defaultHeaderString + ")");
		allParameters.add(OptionBuilder.create(headerShort));
		
		return allParameters;
	}
}



