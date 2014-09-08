package be.svlandeg.diffany.core.listeners;


/**
 * This class implements the {@link ExecutionProgress} interface and prints the progress to 
 * the standard output stream.
 * 
 * @author Sofie Van Landeghem
 */
public class StandardProgressListener implements ExecutionProgress
{
	
	@Override
	public void setProgress(String message, int progress, int total)
	{
		System.out.print(message + ": processed " + progress + " of " + total + " edges.");
	}
	
}
