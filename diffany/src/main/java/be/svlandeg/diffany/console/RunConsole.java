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
		DiffanyOptions dOpt = new DiffanyOptions();
		CommandLineParser parser = new BasicParser();
		
		new RunConsole().runFromConsole(args, dOpt, parser);
	}
	
	/**
	 * Run the Diffany algorithms from the console input, or print the help statement when required/needed.
	 * 
	 * @param args the input arguments provided on the commandline
	 * @param dOpt the {@link DiffanyOptions} object that defines the options of the program
	 * @param parser the {@link CommandLineParser} object that can parse the input arguments
	 */
	public void runFromConsole(String[] args, DiffanyOptions dOpt, CommandLineParser parser)
	{
		Options options = dOpt.getDiffanyOptions();
		
		CommandLine cmd = null;
		try
		{
			cmd = parser.parse(options, args);
		}
		catch(ParseException ex)
		{
			System.err.println( "Could not properly parse the arguments: " + ex.getMessage() );
		}
		
		boolean metaPrinted = displayMetaData(cmd, options);
		
		if (! metaPrinted)
		{
			new RunProject().runAnalysis(cmd);
		}
		
	}
	
	/**
	 * Check whether meta data needs to be displayed, such as the help or version message.
	 * If no request for meta data was needed, this method will return false, signaling the calling method 
	 * that the options object should containc valid project data.
	 * 
	 * @param cmd the parsed input arguments provided on the commandline
	 * @options the options available for running the program
	 * @return whether or not a meta message was printed.
	 */
	private boolean displayMetaData(CommandLine cmd, Options options)
	{
		// TODO: keep version number updated
		String version = "0.0.1";
		
		boolean metaPrinted = false;
		if (cmd.hasOption(DiffanyOptions.helpShort))
		{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("java -jar Diffany_" + version + ".jar", options, true);
			metaPrinted = true;
		}
		else if (cmd.hasOption(DiffanyOptions.versionShort))
		{
			System.out.println("Diffany version: " + version);
			metaPrinted = true;
		}
		return metaPrinted;
	}
	

}