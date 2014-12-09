package be.svlandeg.diffany.examples;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import be.svlandeg.diffany.core.io.NetworkIO;
import be.svlandeg.diffany.core.networks.ConditionNetwork;
import be.svlandeg.diffany.core.networks.ReferenceNetwork;
import be.svlandeg.diffany.core.project.Project;
import be.svlandeg.diffany.core.semantics.DefaultEdgeOntology;
import be.svlandeg.diffany.core.semantics.TreeEdgeOntology;

public class OsmoticUseCase extends GenericExample{


	private final String JAR_DIR = "/data/examples/osmotic_study/input/";
	
	public OsmoticUseCase() 
	{
	}
	
	public Project getTestProject(){
		String name = "Osmotic stress";
		TreeEdgeOntology eo = new DefaultEdgeOntology();
		Project p = new Project(name, eo);
		return p;
	}
	
	public int getTestConfiguration(Project p) throws IOException{
		this.updateStatusMessage("Loading reference network");
		ReferenceNetwork r = getTestReference(JAR_DIR);
		this.updateStatusMessage("Loading condition networks");
		Set<ConditionNetwork> c = getTestConditions(JAR_DIR);
		boolean cleanInput = false;
		this.updateStatusMessage("Adding run configuration to example project");
		int ID = p.addRunConfiguration(r,  c, cleanInput, null);
		this.updateStatusMessage("Example loaded.");
		return ID;
	}

	private ReferenceNetwork getTestReference(String jarDir) throws IOException{
		return NetworkIO.readReferenceNetworkFromResource(jarDir+"Reference network", true);
	}
	
	private Set<ConditionNetwork> getTestConditions(String jarDir) throws IOException{
		Set<ConditionNetwork> cnetworks = new HashSet<ConditionNetwork>();
		this.updateStatusMessage("Loading condition network 1");
		cnetworks.add(getFirstCondition(jarDir));
		this.updateStatusMessage("Loading condition network 2");
		cnetworks.add(getSecondCondition(jarDir));
		this.updateStatusMessage("Loading condition network 3");
		cnetworks.add(getThirdCondition(jarDir));
		this.updateStatusMessage("Loading condition network 4");
		cnetworks.add(getFourthCondition(jarDir));
		
		return cnetworks;
	}

	private ConditionNetwork getFirstCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_1.5h", true);
	}
	private ConditionNetwork getSecondCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_3h", true);
	}
	private ConditionNetwork getThirdCondition(String jarDir) throws IOException{
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_12h", true);
	}	
	private ConditionNetwork getFourthCondition(String jarDir)throws IOException {
		return NetworkIO.readConditionNetworkFromResource(jarDir+"Network_24h", true);
	}

	@Override
	public Project getDefaultProject() {	
		return this.getTestProject();
	}

	@Override
	public int getDefaultRunConfigurationID(Project p) {	
		try {
			return this.getTestConfiguration(p);
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
	}


	private void updateStatusMessage(String message){
		if (this.taskMonitor != null){
			this.taskMonitor.setStatusMessage(message);
		}
	}


	
}
