package be.svlandeg.diffany.cytoscape.gui;

import java.util.Observable;
import java.util.Observer;

import javax.swing.table.AbstractTableModel;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

/**
 * The model for the list of available networks in the selected collection.
 * 
 * @author thpar
 *
 */
public class SelectionTableModel extends AbstractTableModel implements Observer{

	private static final long serialVersionUID = 1L;

	private GUIModel guiModel;

	/**
	 * Column headers
	 */
	String[] columns = {"Include", "Network", "Reference"};
	
	/**
	 * Create a new model based on the general {@link Model} and add it as an {@link Observer} the contained
	 * {@link GUIModel}.
	 * @param model
	 */
	public SelectionTableModel(Model model) {
		this.guiModel = model.getGuiModel();
		this.guiModel.addObserver(this);
	}
	
	
	@Override
	public int getRowCount() {
		return guiModel.getNetworkEntries().size();
	}

	@Override
	public int getColumnCount() {
		return columns.length;
	}

	@Override
	public String getColumnName(int columnIndex) {
		return columns[columnIndex];
	}

	@Override
	public Class<?> getColumnClass(int columnIndex) {
		switch(columnIndex){
		case 0:
		case 2:
			return Boolean.class;
		case 1: 
		default:
			return String.class;
		}
	}

	@Override
	public boolean isCellEditable(int rowIndex, int columnIndex) {
		NetworkEntry entry = guiModel.getNetworkEntries().get(rowIndex);
		switch(columnIndex){
		case 1:
			return false;
		case 0:
			return !entry.isReference();
		case 2:
			return entry.isSelected() && !entry.isReference();
		}
		return false;
	}

	@Override
	public Object getValueAt(int rowIndex, int columnIndex) {
		NetworkEntry entry = guiModel.getNetworkEntries().get(rowIndex);
		switch(columnIndex){
		case 0:
			return entry.isSelected();
		case 1:
			return entry.getName();
		case 2:
			return entry.isReference();
		}
		return false;
	}

	@Override
	public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
		NetworkEntry entry = guiModel.getNetworkEntries().get(rowIndex);
		boolean check = (Boolean)aValue;
		switch(columnIndex){
		case 0:
			entry.setSelected(check);
			break;
		case 2:
			if (!entry.isReference()){
				guiModel.setReferenceEntry(entry);
				this.fireTableDataChanged();
			} 
			break;
		}
	}


	@Override
	public void update(Observable o, Object arg) {
		this.fireTableDataChanged();
	}


}
