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
		    System.out.println("** Version mismatch - Java files don't match library version.");
		    return;
		}
		
		// create engine 
		Rengine re = new Rengine();
		
		// the engine creates R is a new thread, so we should wait until it's ready
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
			System.out.println("Value for java.library.path: " +  System.getProperty( "java.library.path"));
			return;
		}
		new RBridge().getRengine();
	}
}
