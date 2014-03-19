package be.svlandeg.diffany.r;

import org.rosuda.JRI.Rengine;


/** 
 * This class provides an R bridge to allow processing R commands from within Diffany
 * 
 * @author Sofie Van Landeghem
 */
public class RBridge
{

	/**
	 * TODO documentation
	 */
	public void getRengine()
	{
		if (!Rengine.versionCheck()) 
		{
		    System.out.println("Error deploying the Java-R bridge through JRI: Version mismatch of the Rengine Java files.");
		    return;
		}
		
		// create engine 
		Rengine re = new Rengine();
		
		// the engine creates R in a new thread, so we should wait until it's ready
        if (!re.waitForR()) 
        {
            System.out.println("Cannot load R");
            return;
        }
        System.out.println("Succesfully created R engine!");

	}
	
	/**
	 * Currently this method is used for testing
	 * @param args
	 */
	public static void main(String[] args)
	{
		// We have to catch this error to be able to tell the user what to do!
		try
		{
			System.loadLibrary("JRI");
		}
		catch(UnsatisfiedLinkError e)
		{
			System.out.println("Error deploying the Java-R bridge through JRI: " +  e.getMessage());
			System.out.println("");
			System.out.println("This is probably because the directory containing the correct jri.dll is lacking from the path variable.");
			System.out.println("Current value for java.library.path: " +  System.getProperty( "java.library.path"));
			System.out.println("");
			System.out.println("Execution can not continue untill the system settings are fixed!");
			System.out.println("After fixing this issue, please reboot the system and try again.");
			return;
		}
		
		new RBridge().getRengine();
	}
}
