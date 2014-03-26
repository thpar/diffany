package be.svlandeg.diffany.r;

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
		iniEngine(args, runMainLoop, initialCallbacks);
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
	 * Evaluate an R statement through the Java-R bridge.
	 * @return the REXP object resulting from the execution of the R code.
	 */
	protected REXP evaluate(String statement)
	{
		//System.out.println("trying : " + statement);
		return engine.eval(statement);
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
	 * 
	 * @throws IllegalStateException
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
			System.out.println("Error deploying the Java-R bridge through JRI: " + e.getMessage());
			System.out.println("");
			System.out.println("This is probably because the directory containing the correct jri.dll is lacking from the path variable.");
			System.out.println("Current value for java.library.path: " + System.getProperty("java.library.path"));
			System.out.println("");
			System.out.println("Execution can not continue until the system settings are fixed!");
			System.out.println("After fixing this issue, please reboot the system and try again.");
			throw new IllegalStateException("R is not properly installed and configured!");
		}
	}

}
