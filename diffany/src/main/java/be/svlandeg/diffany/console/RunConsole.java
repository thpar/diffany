package be.svlandeg.diffany.console;

import java.io.IOException;

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
	 * First, the commandline is checked for meta options such as 'help' or 'version'. In case of a hit, the appropriate information is printed. If no meta options are provided, the commandline is parsed again for the Diffany arguments.
	 * 
	 * @param args the input arguments provided on the commandline
	 * @param parser the {@link CommandLineParser} object that can parse the input arguments
	 */
	public void runFromConsole(String[] args, CommandLineParser parser)
	{
		Options metaOpt = new MetaOptions().getMetaOptions();
		CommandLine meta_cmd = parseOptions(args, parser, metaOpt, true);

		Options diffOpt = new DiffanyOptions().getDiffanyOptions();
		boolean metaPrinted = displayMetaData(meta_cmd, metaOpt, diffOpt);
		if (!metaPrinted)
		{
			CommandLine diff_cmd = parseOptions(args, parser, diffOpt, false);
			if (diff_cmd != null)
			{
				try
				{
					new RunProject().runAnalysis(diff_cmd);
				}
				catch (IOException ex)
				{
					printError(ex);
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
	 * @param silent if true, ignore exceptions thrown
	 * 
	 * @return the parsed {@link CommandLine} object
	 */
	private CommandLine parseOptions(String[] args, CommandLineParser parser, Options options, boolean silent)
	{
		CommandLine cmd = null;
		try
		{
			cmd = parser.parse(options, args);
		}
		catch (ParseException ex)
		{
			if (! silent)
			{
				printError(ex);
			}
		}

		return cmd;
	}

	/**
	 * Print a user-friendly error message to the console.
	 */
	private void printError(Exception ex)
	{
		String customError = ".";
		if (ex != null)
		{
			customError = ": " + ex.getMessage();
		}
		System.err.println("Problem running Diffany" + customError);
		System.err.println("Run the program with -h or --help (only) to get the help information.");
		System.err.println("Do not add undefined options.");
	}

	/**
	 * Check whether meta data needs to be displayed, such as the help or version message.
	 * 
	 * If no request for meta data was issued, this method will return false, signaling the calling method that the commandline parameters should containc valid project data.
	 * 
	 * @param cmd the parsed input arguments provided on the commandline
	 * @param metaOptions the {@link MetaOptions} defined for this program (e.g. help, version)
	 * @param diffanyOptions the {@link DiffanyOptions} object that defines the options of the actual Diffany program (this information is needed for printing the help message)
	 * 
	 * @return whether or not a meta message was printed.
	 */
	private boolean displayMetaData(CommandLine cmd, Options metaOptions, Options diffanyOptions)
	{
		// TODO: keep version number updated
		String version = "0.0.1";

		boolean metaPrinted = false;

		if (cmd != null && cmd.getOptions() != null && cmd.getOptions().length > 0)
		{
			if (cmd.hasOption(MetaOptions.helpShort))
			{
				HelpFormatter formatter = new HelpFormatter();
				String header = "";
				int consoleSize = getSensibleLength();
				formatter.printHelp(consoleSize, "java -jar Diffany_" + version + ".jar", header, diffanyOptions, header, true);
				metaPrinted = true;
			}
			else if (cmd.hasOption(MetaOptions.versionShort))
			{
				System.out.println("Diffany version: " + version);
				metaPrinted = true;
			}
		}
		return metaPrinted;
	}

	/**
	 * Get an estimate of a proper printing width, depending on the OS
	 * 
	 * @return the length of the string we should be printing
	 */
	private int getSensibleLength()
	{
		String os = System.getProperty("os.name").toLowerCase();

		//Windows
		if (os.contains("win"))
		{
			// Windows Dos Terminal is of width 80 by default
			return 80;
		}
		return 130;
	}

}
