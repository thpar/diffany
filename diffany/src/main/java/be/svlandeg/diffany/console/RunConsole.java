package be.svlandeg.diffany.console;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * This class allows to run the Diffany algorithms from the commandline, by passing the correct arguments through to {@link RunProject}.
 * 
 * @author Sofie Van Landeghem
 */
public class RunConsole
{

	/**
	 * Main method that allows running the Diffany algorithms through the console.
	 * 
	 * @param args the input arguments provided on the commandline
	 */
	public static void main(String[] args)
	{
		CommandLineParser parser = new BasicParser();
		new RunConsole().runFromConsole(args, parser);
	}
	

	/**
	 * Run the Diffany algorithms from the console input, or print the help/version statement when appropriate.
	 * 
	 * First, the commandline is checked for meta options such as 'help' or 'version'. In case of a hit, the appropriate information is printed.
	 * If no meta options are provided, the commandline is parsed again for the Diffany arguments.
	 * 
	 * @param args the input arguments provided on the commandline
	 * @param parser the {@link CommandLineParser} object that can parse the input arguments
	 */
	public void runFromConsole(String[] args, CommandLineParser parser)
	{
		Options metaOpt = new MetaOptions().getMetaOptions();
		CommandLine meta_cmd = parseOptions(args, parser, metaOpt);
		if (meta_cmd != null)
		{
			Options diffOpt = new DiffanyOptions().getDiffanyOptions();
			boolean metaPrinted = displayMetaData(meta_cmd, metaOpt, diffOpt);
			if (!metaPrinted)
			{
				CommandLine diff_cmd = parseOptions(args, parser, diffOpt);
				if (diff_cmd != null)
				{
					new RunProject().runAnalysis(diff_cmd);
				}
			}
		}
	}

	/**
	 * Parse the commandline arguments or print an error message if the parsing fails.
	 * 
	 * @param args the input arguments provided on the commandline
	 * @param parser the {@link CommandLineParser} object that can parse the input arguments
	 * @param mOpt the {@link MetaOptions} defined for this program (e.g. help, version)
	 * 
	 * @return the parsed {@link CommandLine} object
	 */
	private CommandLine parseOptions(String[] args, CommandLineParser parser, Options options)
	{
		CommandLine cmd = null;
		try
		{
			cmd = parser.parse(options, args);
		} catch (ParseException ex)
		{
			System.err.println("Could not properly parse the arguments: " + ex.getMessage());
		}

		return cmd;
	}


	/**
	 * Check whether meta data needs to be displayed, such as the help or version message.
	 * 
	 * If no request for meta data was issued, this method will return false, signaling the calling method 
	 * that the commandline parameters should containc valid project data.
	 * 
	 * @param cmd the parsed input arguments provided on the commandline
	 * @param metaOptions the {@link MetaOptions} defined for this program (e.g. help, version)
	 * @param diffanyOptions the {@link DiffanyOptions} object that defines the options of the actual Diffany program 
	 * (this information is needed for printing the help message)
	 * 
	 * @return whether or not a meta message was printed.
	 */
	private boolean displayMetaData(CommandLine cmd, Options metaOptions, Options diffanyOptions)
	{
		// TODO: keep version number updated
		String version = "0.0.1";

		boolean metaPrinted = false;

		if (cmd.getOptions() != null && cmd.getOptions().length > 0)
		{
			if (cmd.hasOption(MetaOptions.helpShort))
			{
				HelpFormatter formatter = new HelpFormatter();
				//formatter.printHelp("java -jar Diffany_" + version + ".jar", metaOptions, true);
				formatter.printHelp("java -jar Diffany_" + version + ".jar", diffanyOptions, true);
				metaPrinted = true;
			} else if (cmd.hasOption(MetaOptions.versionShort))
			{
				System.out.println("Diffany version: " + version);
				metaPrinted = true;
			}
		}
		return metaPrinted;
	}

}
