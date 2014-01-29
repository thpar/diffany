package be.svlandeg.diffany.console;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
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
		Options options = dOpt.getDiffanyOptions();
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
		try
		{
			cmd = parser.parse(options, args);
		}
		catch(ParseException ex)
		{
			System.err.println( "Could not properly parse the arguments: " + ex.getMessage() );
		}
		
		// TODO catch help/version enquiry here, without sending forward to RunProject
		
		new RunProject().runAnalysis(cmd);
		
	}
	

}