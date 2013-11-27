package be.svlandeg.diffany.cytoscape.gui;

import java.util.Observable;
import java.util.Observer;

import javax.swing.table.AbstractTableModel;

import be.svlandeg.diffany.cytoscape.Model;
import be.svlandeg.diffany.cytoscape.NetworkEntry;

public class SelectionTableModel extends AbstractTableModel implements Observer{

	private static final long serialVersionUID = 1L;

	private GUIModel guiModel;

	String[] columns = {"Include", "Network", "Reference"};
	
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
		switch(columnIndex){
		case 1:
			return false;
		case 0:
			return true;
		case 2:
			NetworkEntry entry = guiModel.getNetworkEntries().get(rowIndex);
			return entry.isSelected();
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
			entry.setReference(check);
			break;
		}
	}


	@Override
	public void update(Observable o, Object arg) {
		this.fireTableDataChanged();
	}


}
