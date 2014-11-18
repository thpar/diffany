package be.svlandeg.diffany.core.progress;

import java.util.Date;


/**
 * This class implements the {@link ExecutionProgress} interface and prints the progress to 
 * the standard output stream, together with a time stamp.
 * 
 * @author Sofie Van Landeghem
 */
public class StandardProgressListener extends ExecutionProgress
{
	
	@Override
	protected void setProgress(String message, int progress, int total)
	{
		System.out.println(new Date() + ": " + message + ": processed " + progress + " of " + total);
	}
	
}