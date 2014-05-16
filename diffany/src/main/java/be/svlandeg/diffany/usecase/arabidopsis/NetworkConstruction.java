package be.svlandeg.diffany.usecase.arabidopsis;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.List;
import java.util.SortedSet;

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
	public void getSignificantGenes(List<OverexpressionData> datasets) throws IOException, URISyntaxException
	{
		GenePrinter gp = new GenePrinter();
		for (OverexpressionData data : datasets)
		{
			System.out.println("");
			System.out.println(data.getName() + ": " + data.getArrayIDs().size() + " IDs");
			SortedSet<String> ids = data.getArrayIDs();
			for (String id : ids)
			{
				double FDR = data.getFDR(id);
				if (FDR <= 0.05)
				{
					if (data.indexedByRawArrayIDs())
					{
						System.out.println("  " + id + " - FDR:" + FDR + " - " + gp.getSynonymsByArrayID(id));
					}
					else
					{
						System.out.println("  " + id + " - FDR:" + FDR + " - " + gp.getSynonymsByLocusID(id));
					}
				}
			}
		}
	}

}
