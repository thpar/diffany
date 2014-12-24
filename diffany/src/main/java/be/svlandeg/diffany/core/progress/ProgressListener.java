package be.svlandeg.diffany.core.progress;

/*
 * #%L
 * Diffany
 * %%
 * Copyright (C) 2014 PSB/UGent - Sofie Van Landeghem and Thomas Van Parys
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU General Lesser Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/lgpl-3.0.html>.
 * #L%
 */


import be.svlandeg.diffany.core.algorithms.CalculateDiff;

/**
 * This interface allows to follow the progress of a differential network calculation,
 * and as such allows to provide feedback to a waiting user.
 * 
 * @author Sofie Van Landeghem
 */
public abstract class ProgressListener
{
	
	private int totalTicks;
	private int doneTicks;

	/**
	 * Iniate a progress listener at a total of 0 ticks.
	 */
	public ProgressListener()
	{
		doneTicks = 0;
		totalTicks = 0;
	}
	
	/**
	 * This method will be called by the differential algorithms in {@link CalculateDiff} every x clicks
	 * to allow a listener to implement e.g. a progress bar or other feedback to the user about the 
	 * duration of the differential algorithm.
	 * 
	 * @param message the message describes the current job that is being done - one ProgressListener can be used for multiple tasks!
	 * @param progress the number of succesful ticks already conducted
	 * @param total the total number of ticks that are envisioned to be needed for full execution
	 */
	protected abstract void setProgress(String message, int progress, int total);
	
	/**
	 * Reset the progress listener to start recording a new task.
	 * @param totalTicks the total number of ticks the new task will take
	 */
	public void reset(int totalTicks)
	{
		doneTicks = 0;
		this.totalTicks = totalTicks;
	}
	
	/**
	 * Report the fact that there has been some progress, expressed by a number of ticks that were performed
	 * 
	 * @param message the message that should be displayed to the user at this point, i.e. the task that is now being done
	 * @param ticks the number of ticks that have passed since the last 'progress update'
	 * 
	 */
	public void addTicks(String message, int ticks)
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
		setProgress(message, doneTicks, totalTicks);
	}
	
	/**
	 * Obtain the number of ticks that are still to be done for this task to complete
	 * @return the number of ticks that still need to be done
	 */
	protected int ticksToGo()
	{
		return totalTicks - doneTicks;
	}
	
}
