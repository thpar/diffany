package be.svlandeg.diffany.core.listeners;

import java.util.Date;


/**
 * This class implements the {@link ExecutionProgress} interface and prints the progress to 
 * the standard output stream, together with a time stamp.
 * 
 * @author Sofie Van Landeghem
 */
public class StandardProgressListener implements ExecutionProgress
{
	
	@Override
	public void setProgress(String message, int progress, int total)
	{
		System.out.println(new Date() + ": " + message + ": processed " + progress + " of " + total);
	}
	
}
