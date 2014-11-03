package be.svlandeg.diffany.examples;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.DefaultNodeMapper;
import be.svlandeg.diffany.core.semantics.NodeMapper;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

public class OsmoticUseCase extends GenericExample{

	private NodeMapper nm;

	private final String JAR_DIR = "/data/osmotic/";
	
	public OsmoticUseCase() {
		nm = new DefaultNodeMapper();
	}
	
	public Project getTestProject(){
		String name = "Osmotic Use Case";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo, nm);
		return p;
	}
	
	public int getTestConfiguration(Project p) throws IOException{
		ReferenceNetwork r = getTestReference(JAR_DIR);
		Set<ConditionNetwork> c = getTestConditions(JAR_DIR);
		boolean cleanInput = true;
		int ID = p.addRunConfiguration(r,  c, cleanInput, null);
		return ID;
	}

	private ReferenceNetwork getTestReference(String jarDir) throws IOException{
		return NetworkIO.readReferenceNetworkFromResource(jarDir+"Reference network", nm, true);
	}
	
	private Set<ConditionNetwork> getTestConditions(String jarDir) throws IOException{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();
		cnetworks.add(getFirstCondition(jarDir));
		cnetworks.add(getSecondCondition(jarDir));
		cnetworks.add(getThirdCondition(jarDir));
		cnetworks.add(getFourthCondition(jarDir));
		return cnetworks;
	}

	private ConditionNetwork getFirstCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_1.5h", nm, true);
	}
	private ConditionNetwork getSecondCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_3h", nm, true);
	}
	private ConditionNetwork getThirdCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_12h", nm, true);
	}	
	private ConditionNetwork getFourthCondition(String jarDir)throws IOException {
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_24h", nm, true);
	}




	
}
