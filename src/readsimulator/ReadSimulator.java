package readsimulator;

import java.util.HashMap;

public class ReadSimulator implements Runnable {

	private HashMap<String, HashMap<String, Integer>> readcounts;
	private HashMap<String, int[]> fasta_annotation;

	public ReadSimulator() {

		readcounts = new HashMap<>();


	}

	@Override
	public void run() {

		readFiles();

		printHeaders();

		calculateFragments();


	}

	private void readFiles() {
		/*
		 * Already done by FileReader: get Paths and inputParameters
		 * TODO: Read Readcounter.simulation
		 * TODO: Read GTF File, but only lines with "CDS" and ID's from Readcounter.simulation
		 * TODO: Read FASTA_indexfile
		 * */




	}

	private void printHeaders() {

	}

	private void calculateFragments() {

	}

}

