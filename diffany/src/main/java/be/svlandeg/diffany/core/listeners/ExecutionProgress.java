package be.svlandeg.diffany.core.listeners;

import be.svlandeg.diffany.core.algorithms.CalculateDiff;

/**
 * This interface allows to follow the progress of a differential network calculation,
 * and as such allows to provide feedback to a waiting user.
 * 
 * @author Sofie Van Landeghem
 */
public interface ExecutionProgress
{
	
	/**
	 * This method will be called by the differential algorithms in {@link CalculateDiff} every x clicks
	 * to allow a listener to implement e.g. a progress bar or other feedback to the user about the 
	 * duration of the differential algorithm.
	 * 
	 * @param message the message describes the current job that is being done - one ExecutionProgress can be used for multiple tasks!
	 * @param progress the number of succesful ticks already conducted
	 * @param total the total number of ticks that are envisioned to be needed for full execution
	 */
	public void setProgress(String message, int progress, int total);
	
}
