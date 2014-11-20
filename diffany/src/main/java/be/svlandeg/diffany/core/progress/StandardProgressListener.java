package be.svlandeg.diffany.core.progress;

import java.util.Date;


/**
 * This class implements the {@link ProgressListener} interface and prints the progress to 
 * the standard output stream, together with a time stamp.
 * 
 * @author Sofie Van Landeghem
 */
public class StandardProgressListener extends ProgressListener
{
	
	private boolean print;
	
	/**
	 * Constructor which specifies whether or not to print the progress messages.
	 * Usually, this value will be true, except for the JUnit usage of this class.
	 * @param print whether or not to print the progress messages
	 */
	public StandardProgressListener(boolean print)
	{
		super();
		this.print = print;
	}
	
	@Override
	protected void setProgress(String message, int progress, int total)
	{
		if (print)
		{
			System.out.println(new Date() + ": " + message + ": processed " + progress + " of " + total);
		}
	}
	
}
