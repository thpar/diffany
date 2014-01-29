package be.svlandeg.diffany.console;

import org.apache.commons.cli.CommandLine;

/**
 * This class can run the Diffany algorithms from a {@link org.apache.commons.cli.CommandLine} object.
 * 
 * @author Sofie Van Landeghem
 */
public class RunProject
{

	/**
	 * Run a Diffany analysis, depending on the parameters provided on the commandline.
	 * 
	 * @param cmd the CommandLine object containing all parameters provided by the user
	 * @throws IllegalArgumentException when a crucial argument is missing
	 */
	public void runAnalysis(CommandLine cmd)
	{
		boolean toLog = false;
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			toLog = true;
		}
		
		String refDir = cmd.getOptionValue(DiffanyOptions.refShort);
		if (refDir == null)
		{
			throw new IllegalArgumentException("Fatal error: please provide a valid refDir argument!");
		}
		
		String condDir = cmd.getOptionValue(DiffanyOptions.conShort);
		if (condDir == null)
		{
			throw new IllegalArgumentException("Fatal error: please provide a valid refDir argument!");
		}
		
		
		
		// TEST CODE
		if (toLog)
		{
			System.out.println("Will log");
		}
		else
		{
			System.out.println("Will not log");
		}

	}

}
