package be.svlandeg.diffany.core.progress;

/**
 * A scheduled task knows how many ticks it will take to execute in total. This value should be used to prompt the ProgressListener at regular intervals about the progress.
 * 
 * @author Sofie Van Landeghem
 */
public class ScheduledTask
{
	
	private int totalTicks;
	private String message;
	private ExecutionProgress ep;
	private int doneTicks;

	/**
	 * Constructor
	 * @param ep TODO
	 * @param totalTicks TODO
	 */
	public ScheduledTask(ExecutionProgress ep, int totalTicks)
	{
		this.ep = ep;
		this.totalTicks = totalTicks;
		doneTicks = 0;
		message = "Calculating ...";
	}
	
	public void setMessage(String message)
	{
		this.message = message;
	}
	
	public void ticksDone(int ticks)
	{
		if (ticks > ticksToGo())
		{
			throw new IllegalArgumentException("Can not report more progress than what was scheduled!");
		}
		if (ticks < 1)
		{
			throw new IllegalArgumentException("Can not report " + ticks + " ticks!");
		}
		doneTicks += ticks;
		ep.addTicks(message, ticks);
	}
	
	public void done()
	{
		if (ticksToGo() > 0)
		{
			ep.addTicks(message, ticksToGo());
			doneTicks = totalTicks;
		}
	}
	
	public int ticksToGo()
	{
		return totalTicks - doneTicks;
	}
	
	protected int totalTicks()
	{
		return totalTicks;
	}
	
	public ExecutionProgress getListener()
	{
		return ep;
	}
	
	
}
