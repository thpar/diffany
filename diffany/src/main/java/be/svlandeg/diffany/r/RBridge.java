package be.svlandeg.diffany.r;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

/**
 * This class provides an R bridge to allow processing R commands from within Diffany.
 * The R commands can be defined in-line, or in a separate .R file.
 * Once the RBridge is created, use the ExecuteR class to actually execute statements!
 * 
 * @author Sofie Van Landeghem
 */
public class RBridge
{

	private final static String[] DEFAULT_ARGS = new String[]{"--vanilla"};
	private final static boolean DEFAULT_LOOP = false;
	private Rengine engine;

	private String logfile;

	/**
	 * Produces a bridge to R. The bridge should be closed when the application is done with the R engine!
	 * 
	 * @param args arguments to be passed to R (not an empty list!)
	 * @param runMainLoop if set to true the the event loop will be started as soon as possible, otherwise no event loop is started
	 * @param initialCallbacks an instance implementing the RMainLoopCallbacks interface that provides methods to be called by R (can be null)
	 * @param logfile the absolute file path which will contain the R log file (if null, will be ignored). Will be emptied first.
	 */
	public RBridge(String[] args, boolean runMainLoop, RMainLoopCallbacks initialCallbacks, String logfile)
	{
		checkSystem();
		iniEngine(args, runMainLoop, initialCallbacks);
		if (logfile != null)
		{
			this.logfile = logfile;
			try
			{
				// Make sure the file is empty before starting to write the log
				FileWriter writer = new FileWriter(logfile, false);
				writer.close();
				divertLog(logfile);
			}
			catch (IOException ex)
			{
				System.out.println(" ! error: could not divert the output to " + logfile);
				endLogDiversion();
			}
		}
	}

	/**
	 * Produces a bridge to R. The bridge should be closed when the application is done with the R engine!
	 * This constructor will call the main constructor with vanilla arguments, runMainLoop=false and no initialCallBacks (null).
	 * 
	 * @param logfile the absolute file path which will contain the R log file
	 */
	public RBridge(String logfile)
	{
		this(DEFAULT_ARGS, DEFAULT_LOOP, null, logfile);
	}

	/**
	 * Produces an engine for R, and initializes this engine as a protected field within this class.
	 * The function end should be called when the application is done with the engine!
	 * 
	 * @param args arguments to be passed to R (not an empty list!)
	 * @param runMainLoop if set to true the the event loop will be started as soon as possible, otherwise no event loop is started
	 * @param initialCallbacks an instance implementing the RMainLoopCallbacks interface that provides methods to be called by R (can be null)
	 * 
	 * @throws IllegalStateException when something went wrong while creating the R engine, usually due to misconfiguration.
	 */
	private void iniEngine(String[] args, boolean runMainLoop, RMainLoopCallbacks initialCallbacks) throws IllegalStateException
	{
		if (!Rengine.versionCheck())
		{
			System.out.println("Error deploying the Java-R bridge through JRI: Version mismatch of the Rengine Java files.");
			throw new IllegalStateException("R is not properly installed and configured!");
		}

		engine = new Rengine(args, runMainLoop, initialCallbacks);

		if (!engine.waitForR())
		{
			System.out.println("Error deploying the Java-R bridge through JRI: cannot load R");
			throw new IllegalStateException("R is not properly installed and configured!");
		}
	}

	/**
	 * R wants only forward slashes in a path location, so this method converts backward slashes to forward.
	 * 
	 * @param dir_path the original path, which may contain both forward and backward slashes
	 * @return a path location usable by the R bridge (i.e. containing only forward slashes)
	 */
	public String convertSlashes(String dir_path)
	{
		return dir_path.replace("\\", "/");
	}

	/**
	 * Evaluate an R statement through the Java-R bridge.
	 * @param statement the original R statement
	 * 
	 * @return the REXP object resulting from the execution of the R code.
	 */
	protected REXP evaluate(String statement)
	{
		//System.out.println("trying : " + statement);
		return engine.eval(statement);
	}

	/**
	 * Print the logs of the R environment into the specified file. The file is appended and NOT cleaned beforehand!
	 * 
	 * @param filepath the absolute file path which will contain the R log file
	 */
	private void divertLog(String filepath)
	{
		evaluate("log<-file('" + convertSlashes(filepath) + "')");
		evaluate("sink(log, append=TRUE)");
		evaluate("sink(log, append=TRUE, type='message')");
	}

	/**
	 * Stop printing the R outputs to the file specified earlier.
	 * This method will also close all existing connections to ensure proper cleanup.
	 */
	private void endLogDiversion()
	{
		evaluate("sink(file = NULL)");
		evaluate("closeAllConnections()");
	}

	/**
	 * Close the R engine: should be called when the application is done with R calculations. Calls an asynchronous method!
	 * This method will also close the connection to the log file, if any.
	 */
	public void close()
	{
		engine.end();
		endLogDiversion();
	}

	/**
	 * Retrieve all errors found in the logfile of the R commando's. Specifically, this will return all lines starting with "Error".
	 * When there was a problem reading the logfile, this will be specified as an error message in the returned list itself.
	 * 
	 * @return a chronological list of errors which happened during execution of the R commando's through this bridge.
	 */
	public List<String> getErrorsFromLogfile()
	{
		List<String> errors = new ArrayList<String>();
		if (logfile == null)
		{
			errors.add("No logfile specified!");
		}
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(logfile));
			String line = reader.readLine();
			while (line != null)
			{
				if (line.startsWith("Error"))
				{
					errors.add(line);
				}
				line = reader.readLine();
			}
			reader.close();
		}
		catch(IOException ex)
		{
			errors.add("Error reading logfile: " + ex.getMessage());
		}
		
		return errors;
	}

	/**
	 * This private message will inform the user of inappropriately installed R/JRI versions,
	 * and provide help information and pointers for fixing the issue.
	 * 
	 * @throws IllegalStateException when R is not properly installed or configured
	 */
	private void checkSystem() throws IllegalStateException
	{
		/* We have to catch this error to be able to tell the user what to do! */
		try
		{
			System.loadLibrary("JRI");
		}
		catch (UnsatisfiedLinkError e)
		{
			System.err.println(" ! Error deploying the Java-R bridge through JRI: " + e.getMessage());
			System.err.println("");

			System.err.println(" ! This may be because the R installation is not on the path variable (e.g. C:/Program Files/R/R-X.Y.Z/bin/x64/)");
			System.err.println(" ! Or because the jri.dll is not on the path variable (e.g. C:/Program Files/R/R-X.Y.Z/library/rJava/jri/x64/)");
			System.err.println("");
			System.err.println(" ! If R is not installed, download and install it from http://www.r-project.org/");
			System.err.println(" ! If rJava/JRI is not installed, run this command in R: << install.packages(\"rJava\") >>");
			System.err.println(" ! A nice GUI for R can be found at http://www.rstudio.com/");
			System.err.println("");
			System.err.println(" ! Current value for java.library.path: " + System.getProperty("java.library.path"));
			System.err.println(" ! This can be edited by changing the system environment variables.");
			System.err.println("");
			System.err.println(" ! Execution can not continue until the system settings are fixed!");
			System.err.println(" ! After fixing this issue, please reboot the system and try again.");
			System.err.println("");
			throw new IllegalStateException("R is not properly installed and configured!");
		}
	}

}
