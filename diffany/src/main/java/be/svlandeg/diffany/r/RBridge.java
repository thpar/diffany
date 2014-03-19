package be.svlandeg.diffany.r;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RMainLoopCallbacks;
import org.rosuda.JRI.Rengine;

/** 
 * This class provides an R bridge to allow processing R commands from within Diffany.
 * The R commands can be defined in-line, or in a separate .R file.
 * 
 * @author Sofie Van Landeghem
 */
public class RBridge
{

	protected final static String[] DEFAULT_ARGS = new String[]{ "--vanilla" };
	protected final static boolean DEFAULT_LOOP = false;
	protected Rengine engine;

	/**
	 * Produces a bridge to R. The bridge should be closed when the application is done with the R engine!
	 * 
	 * @param args arguments to be passed to R (not an empty list!)
	 * @param runMainLoop if set to true the the event loop will be started as soon as possible, otherwise no event loop is started
	 * @param initialCallbacks an instance implementing the RMainLoopCallbacks interface that provides methods to be called by R (can be null)
	 */
	public RBridge(String[] args, boolean runMainLoop, RMainLoopCallbacks initialCallbacks)
	{
		checkSystem();
		engine = getRengine(args, runMainLoop, initialCallbacks);
	}

	/**
	 * Produces a bridge to R. The bridge should be closed when the application is done with the R engine!
	 * This constructor will call the main constructor with vanilla arguments, runMainLoop=false and no initialCallBacks (null).
	 * 
	 * @param args arguments to be passed to R (not an empty list!)
	 * @param runMainLoop if set to true the the event loop will be started as soon as possible, otherwise no event loop is started
	 * @param initialCallbacks an instance implementing the RMainLoopCallbacks interface that provides methods to be called by R (can be null)
	 */
	public RBridge()
	{
		this(DEFAULT_ARGS, DEFAULT_LOOP, null);
	}

	/**
	 * Produces an engine for R. The function end should be called when the application is done with the engine.
	 * 
	 * @param args arguments to be passed to R (not an empty list!)
	 * @param runMainLoop if set to true the the event loop will be started as soon as possible, otherwise no event loop is started
	 * @param initialCallbacks an instance implementing the RMainLoopCallbacks interface that provides methods to be called by R (can be null)
	 * 
	 * @return the R engine (should not be null)
	 * @throws IllegalStateException when something went wrong while creating the R engine, usually due to misconfiguration.
	 */
	protected Rengine getRengine(String[] args, boolean runMainLoop, RMainLoopCallbacks initialCallbacks) throws IllegalStateException
	{
		if (!Rengine.versionCheck())
		{
			System.out.println("Error deploying the Java-R bridge through JRI: Version mismatch of the Rengine Java files.");
			throw new IllegalStateException("R is not properly installed and configured!");
		}

		Rengine re = new Rengine(args, runMainLoop, initialCallbacks);

		if (!re.waitForR())
		{
			System.out.println("Error deploying the Java-R bridge through JRI: cannot load R");
			throw new IllegalStateException("R is not properly installed and configured!");
		}
		return re;
	}

	/**
	 * Testing method : R instructions from file
	 * TODO: will this code work when packaged inside a jar or will we need to create a tmp file?
	 */
	public void randomTesting1() 
	{
		System.out.println("Test: reading R script from file");
		
		URL scriptURL = Thread.currentThread().getContextClassLoader().getResource("helloWorld.R");

		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(new File(scriptURL.toURI())));
			String line = reader.readLine();
			while (line != null)
			{
				engine.eval(line);
				line = reader.readLine();
			}
			reader.close();
		}
		catch (IOException e)
		{
			System.out.println("Couldn't read R code : " + e.getMessage());
			return;
		}
		catch (URISyntaxException e)
		{
			System.out.println("Couldn't read R code : " + e.getMessage());
			return;
		}

		REXP string_greeting = engine.eval("greeting");
		String s = string_greeting.asString();
		
		REXP number_greeting = engine.eval("number");
		double i = number_greeting.asDouble();
		
		System.out.println(" Greeting from R: " + s + " - " + i);
	}

	/**
	 * Testing method: in-code R
	 */
	public void randomTesting2()
	{
		System.out.println("Test: executing R code directly in Java code");
		
		engine.eval(String.format("greeting <- '%s'", "Hello R In-code World"));
		REXP result = engine.eval("greeting");
		System.out.println(" Greeting from R: " + result.asString());
	}

	/**
	 * Close the R engine, should be called when the application is done with R calculations.
	 * Calls an asynchronous method!
	 */
	public void close()
	{
		engine.end();
	}

	/**
	 * This private message will inform the user of inappropriately installed R/JRI versions.
	 * @throws IllegalStateException
	 */
	private void checkSystem() throws IllegalStateException
	{
		// We have to catch this error to be able to tell the user what to do!
		try
		{
			System.loadLibrary("JRI");
		}
		catch (UnsatisfiedLinkError e)
		{
			System.out.println("Error deploying the Java-R bridge through JRI: " + e.getMessage());
			System.out.println("");
			System.out.println("This is probably because the directory containing the correct jri.dll is lacking from the path variable.");
			System.out.println("Current value for java.library.path: " + System.getProperty("java.library.path"));
			System.out.println("");
			System.out.println("Execution can not continue untill the system settings are fixed!");
			System.out.println("After fixing this issue, please reboot the system and try again.");
			throw new IllegalStateException("R is not properly installed and configured!");
		}
	}

	/**
	 * Currently this method is used for testing
	 */
	public static void main(String[] args) 
	{
		RBridge bridge = new RBridge();
		
		bridge.randomTesting1();
		System.out.println(" ");
		
		bridge.randomTesting2();
		System.out.println(" ");
		
		bridge.close();
	}
}
