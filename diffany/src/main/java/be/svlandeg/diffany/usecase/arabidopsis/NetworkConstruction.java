package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class allows to construct networks out of overexpression/coexpression values.
 * 
 * @author Sofie Van Landeghem
 */
public class NetworkConstruction
{

	/**
	 * TODO
	 */
	public NetworkConstruction()
	{

	}

	/**
	 * 
	 * @param datasets
	 * @throws URISyntaxException 
	 * @throws IOException 
	 */
	public void getSignificantGenes(List<OverexpressionData> datasets, double threshold) throws IOException, URISyntaxException
	{
		for (OverexpressionData data : datasets)
		{
			System.out.println("");
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs analysed");
			SortedSet<String> ids = data.getArrayIDs();
			SortedSet<String> sign_ids_up = new TreeSet<String>();
			SortedSet<String> sign_ids_down = new TreeSet<String>();
			for (String id : ids)
			{
				double FDR = data.getFDR(id);
				if (FDR <= threshold)
				{
					double FC = data.getFoldchange(id);
					if (FC > 0)
					{
						sign_ids_up.add(id);
					}
					else
					{
						sign_ids_down.add(id);
					}
				}
			}
			System.out.print("  " + sign_ids_up.size() + " IDs upregulated: ");
			for (String id : sign_ids_up)
			{
				System.out.print("  " + id);
			}
			System.out.println("");
			System.out.print("  " + sign_ids_down.size() + " IDs downregulated: ");
			for (String id : sign_ids_down)
			{
				System.out.print("  " + id);
			}
			System.out.println("");
		}
	}

}
