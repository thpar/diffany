package be.svlandeg.diffany.console;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;



/**
 * This class allows to run the Diffany algorithms from the commandline, by parsing the correct arguments through to {@link RunProject}.
 * 
 * @author Sofie Van Landeghem
 */
public class RunConsole
{
	
	/**
	 * Main method that allows running the Diffany algorithms through the console.
	 * 
	 * @param args
	 * @throws ParseException when the arguments can not be parsed properly (TODO friendly error)
	 */
	public static void main(String[] args) throws ParseException 
	{
		CommandLineParser parser = new BasicParser();
		DiffanyOptions dOpt = new DiffanyOptions();
		Options options = dOpt.getDiffanyOptions();
		
		CommandLine cmd = parser.parse(options, args);
		new RunProject().runAnalysis(cmd);
	}
	
	

}
