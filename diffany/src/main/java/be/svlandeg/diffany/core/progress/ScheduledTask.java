package be.svlandeg.diffany.core.progress;

/**
 * A scheduled task knows how many ticks it will take to execute in total. 
 * This value should be used to prompt the ProgressListener at regular intervals about the progress.
 * 
 * @author Sofie Van Landeghem
 */
public class ScheduledTask
{
	
	private int totalTicks;
	private String message;
	private ProgressListener ep;
	private int doneTicks;

	/**
	 * Constructor: create a new scheduled task belonging to a certain progress listener and with a pre-defined number of total ticks.
	 * 
	 * @param ep the progress listener
	 * @param totalTicks the total number of steps required for this task
	 */
	public ScheduledTask(ProgressListener ep, int totalTicks)
	{
		this.ep = ep;
		this.totalTicks = totalTicks;
		doneTicks = 0;
		message = "Calculating ...";
	}
	
	/**
	 * Set the message belonging to this task
	 * @param message the message that should tell the user what this task does
	 */
	public void setMessage(String message)
	{
		this.message = message;
	}
	
	/**
	 * Report that a number of ticks has been done since the last report.
	 * @param ticks the number of ticks that were done
	 */
	public void ticksDone(int ticks)
	{
		if (ticks > ticksToGo())
		{
			ticks = ticksToGo();
		}
		if (ticks > 0)
		{
			doneTicks += ticks;
			ep.addTicks(message, ticks);
		}
	}
	
	/**
	 * Report that the task is finished.
	 */
	public void done()
	{
		if (ticksToGo() > 0)
		{
			ep.addTicks(message, ticksToGo());
			doneTicks = totalTicks;
		}
	}
	
	/**
	 * Return the number of ticks that still need to be done to succesfully complete this task.
	 * @return the number of ticks that still need to be done
	 */
	public int ticksToGo()
	{
		return totalTicks - doneTicks;
	}
	
	/**
	 * Return the total number of original ticks in this task (not taking into account how many of them have already been done)
	 * @return the total number of original ticks in this task
	 */
	protected int totalTicks()
	{
		return totalTicks;
	}
	
	/**
	 * Retrieve the progress listener object
	 * @return the progress listener that reports the progress to the user
	 */
	public ProgressListener getListener()
	{
		return ep;
	}
	
	
}
