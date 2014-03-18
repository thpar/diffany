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
		    System.exit(1);
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
	
	public static void main(String[] args)
	{
		new RBridge().getRengine();
	}
}
