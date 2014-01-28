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
	 * @param cmd
	 */
	public void runAnalysis(CommandLine cmd)
	{
		boolean toLog = false;
		if (cmd.hasOption(DiffanyOptions.logShort))
		{
			toLog = true;
		}
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
