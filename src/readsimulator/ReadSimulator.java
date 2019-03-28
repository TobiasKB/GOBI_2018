package readsimulator;

import exonskipping.Gen;
import exonskipping.Transcript;
import utils.CommandLine_Parser;
import utils.GTF_FileParser;
import utils.TesException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

public class ReadSimulator implements Runnable {

	/*
	 * readcounts HashMap<gene_id < transcript_id , count_for_transcript>
	 * fasta_annotation HashMap < chr, [entry_length, entry_start, line_length, line_length_nlChar]>
	 * target_genes < gene_id, Gen >
	 * */
	private HashMap<String, HashMap<String, Integer>> readcounts;
	private HashMap<String, long[]> fasta_annotation;
	private HashMap<String, Gen> target_genes;

	public ReadSimulator() {

		readcounts = new HashMap<String, HashMap<String, Integer>>();
		fasta_annotation = new HashMap<String, long[]>();
		target_genes = new HashMap<>();

	}

	@Override
	public void run() {

		readFiles();

		printHeaders();

		getSequences();

		calculateFragments();


	}

	private void readFiles() {
		/*
		 * Already done by FileReader: get Paths and inputParameters
		 * Ein Thread liest das IndexFile , einer den Rest//
		 *
		 * */


		new Thread(new Runnable() {
			public void run() {
				try {
//				Auslesen des FastaindexFiles
					BufferedReader br = new BufferedReader(new FileReader(CommandLine_Parser.inputFile_fidx));
					String line = br.readLine();
					while (line != null) {
						String[] splitline = line.split("\t");
						if (!fasta_annotation.containsKey(splitline[0])) {
							fasta_annotation.put(splitline[0], new long[]{Long.parseLong(splitline[1]), Long.parseLong(splitline[2]), Long.parseLong(splitline[3]), Long.parseLong(splitline[4])});
						}

						line = br.readLine();
					}
					br.close();

				} catch (IOException e) {
					throw new TesException("Could not Read inputFile_fidx", e);
				}
			}
		}).start();

		new Thread(new Runnable() {
			public void run() {
				try {
//				Auslesen des FastaindexFiles
					BufferedReader br = new BufferedReader(new FileReader(CommandLine_Parser.inputFile_readcounts));
					String line = br.readLine();
					while (line != null) {
						if (line.contains("gene")) {
							line = br.readLine();
							continue;
						}
						String[] splitline = line.split("\t");

						if (!readcounts.containsKey(splitline[0])) {

							HashMap h = new HashMap<String, Integer>();
							h.put(splitline[1], Integer.valueOf(splitline[2]));

							readcounts.put(splitline[0], h);
						} else {
							readcounts.get(splitline[0]).put(splitline[1], Integer.valueOf(splitline[2]));
						}
						line = br.readLine();
					}
				} catch (IOException e) {
					throw new TesException("Fehler beim einlesen des Fastaindexfiles fidx", e);
				}
			}
		}).start();

		try {

			BufferedReader br = new BufferedReader(new FileReader(CommandLine_Parser.inputFile_GTF));
			String line = br.readLine();
			String[] cleanline;

			while (line != null) {
				cleanline = GTF_FileParser.parseLine(line, null, null);
				if (cleanline != null) {

					for (String gene_id : readcounts.keySet()) {

						for (String transcript_id : readcounts.get(gene_id).keySet()) {

							if (cleanline[5].equals(gene_id) && cleanline[7].equals(transcript_id)) {

								if (!target_genes.containsKey(cleanline[5])) {

									target_genes.put(cleanline[5], new Gen(cleanline[5], cleanline[4], cleanline[0], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[3], cleanline[6]));

								}
								if (!target_genes.get(cleanline[5]).get_Transcripts().containsKey(cleanline[7]))
									target_genes.get(cleanline[5]).add_Transcript(new Transcript(cleanline[7], cleanline[3], cleanline[6], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[5]));

								target_genes.get(cleanline[5]).add_Region(cleanline[7], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]));
								target_genes.get(cleanline[5]).add_nprots();
							}
						}
					}
				}
				line = br.readLine();
			}

			br.close();

		} catch (IOException e) {

		}
	}

	private void getSequences() {
		/*
		 * TODO: Jump into FASTA file and grab Sequences of relevance /Safe Sequence to Transcript
		 * Going through all Transcript IDs for all Genes, with the gene giving information about the position in the genome (chr)
		 * For each perform a lookup in the fidx file (saved in "fasta_annotation" to gather the position more quickly.
		 *
		 * For final access, a random access file is used, using the java library
		 * */
		try {


			RandomAccessFile raf = new RandomAccessFile(CommandLine_Parser.inputFile_fasta, "r");

			for (Map.Entry<String, Gen> gen : target_genes.entrySet()) {

				gen.getValue().get_Transcripts().forEach((transcript_id, t) -> {
					StringBuilder stringBuilder = new StringBuilder();
/*
					System.out.println(t.getTrans_id());
					System.out.println(t.getStart());
					System.out.println(t.getEnd());*/
					char[] temp_sequence = new char[t.getEnd() - t.getStart()];


					/*
					 *TODO: For Schleife ueber alle Exons, nicht von Transcript Start bis Ende und zusammencutten.
					 *Noch nicht ganz korrekt, leichte abweichung, nicht ganz so viele gefunden wie eigentlich benoetigt. Ausserdem falscher Index, somewhere, weil immer noch kein Match .
					 * */

					t.getExons().forEach(exon -> {
						try {
							raf.seek((fasta_annotation.get(gen.getValue().getchr())[1] + exon.getStart()));
							byte[] bytey = new byte[(exon.getEnd() - exon.getStart())];

							raf.readFully(bytey, 0, exon.getEnd() - exon.getStart());

							String temp = new String(bytey, StandardCharsets.UTF_8)/*.replaceAll("\n", "")*/;

							stringBuilder.append(temp);


						} catch (IOException e) {
							throw new TesException("Error Seeking correnct Line in FastaFile", e);
						}
					});

					System.out.println(stringBuilder.toString());
					System.out.println();
				});
			}

		} catch (IOException e) {
			throw new TesException("Failed to read Fasta File @RandomFileAccess", e);
		}

	}

	private void calculateFragments() {


	}

	private void printHeaders() {

//		TODO: Print Headers to outputfiles; generate outputfiles. Remember to concat, not overwrite them later on!


	}


}

