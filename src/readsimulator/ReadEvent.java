package readsimulator;

import java.util.List;

public class ReadEvent {

	private String id;
	private String chr;
	private String gene;
	private String fw_sequence;
	private String rw_sequence;
	private List<Integer> fw_mutations;
	private List<Integer> rw_mutations;
	private int[] t_fw_coordinates;
	private int[] t_rw_coordinates;
	private List<Integer> fw_coordinates;
	private List<Integer> rw_coordinates;

	public ReadEvent() {


	}


	@Override
	public String toString() {
		return " ";
	}

}
