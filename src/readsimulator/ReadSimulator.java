package readsimulator;

import exonskipping.Gen;
import exonskipping.Transcript;
import org.apache.commons.math3.distribution.NormalDistribution;
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
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import static utils.CommandLine_Parser.*;

public class ReadSimulator implements Runnable {

	/*
	 * readcounts HashMap<gene_id < transcript_id , count_for_transcript>
	 * fasta_annotation HashMap < chr, [entry_length(0), entry_start(1), line_length(2), line_length_nlChar(3)]>
	 * target_genes < gene_id, Gen >
	 * */
	private HashMap<String, HashMap<String, Integer>> readcounts;
	private HashMap<String, long[]> fasta_annotation;
	private HashMap<String, Gen> target_genes;
	/*
	 * @fw_reads: < read_id, sequence>
	 * @rw_reads: < read_id, sequence>
	 * */
	private HashMap<Integer, String> fw_reads;
	private HashMap<Integer, String> rw_reads;

	/*
	 * @fw_mutations:<read_id,  mutations>
	 * @rw_mutations:<read_id,  mutations>
	 * */
	private HashMap<Integer, String> fw_mutations_pointer;
	private HashMap<Integer, String> rw_mutations_pointer;

	/*
	 * Liste mit allen Events pro Readid//Mappinginfo File
	 * */
	private HashMap<String, String> read_mappinginfo;

	/*
	 * fw_regvec position des fw relativ auf dem Genom/Chromosom --> Transkript besteht aus mehr als nur Exons;
	 * allerdings werden nur Exons betrachtet --> gib die Position der jeweiligen Exons aus.
	 * Also fur Transkript fw 0-75 kann auf 2 oder mehr Exons aufgeteilt werden.
	 * --> get Exon? get Exon position?
	 * rw_regvec
	 *
	 * t_fw_regvec position des fw relativ auf dem Transkript
	 * t_rw regvec positopm des rw relativ auf dem Transkript
	 * Beispiel: Transkript beginnt bei 390842009; Laenge: 305 // ein Read liegt bei 0-75, einer bei 230-305
	 * (0 based, End EXCLUSIVE)
	 * */
	private NormalDistribution rand;


	public ReadSimulator() {

		readcounts = new HashMap<String, HashMap<String, Integer>>();
		fasta_annotation = new HashMap<String, long[]>();
		target_genes = new HashMap<>();
		fw_mutations_pointer = new HashMap<>();
		rw_mutations_pointer = new HashMap<>();
		read_mappinginfo = new HashMap<>();
		fw_reads = new HashMap<>();
		rw_reads = new HashMap<>();


	}

	@Override
	public void run() {

		readFiles();

		printHeaders();

		getSequences();

		calculate_core();

		print();

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
//				Auslesen des reascount files
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
							if (readcounts.get(splitline[0]).containsKey(splitline[1])) {
								readcounts.get(splitline[0]).put(splitline[1], readcounts.get(splitline[0]).get(splitline[1] + Integer.valueOf(splitline[2])));
							} else {
								readcounts.get(splitline[0]).put(splitline[1], Integer.valueOf(splitline[2]));
							}
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

								target_genes.get(cleanline[5]).add_Region(cleanline[7], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[8]);
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
		 * Going through all Transcript IDs for all Genes, with the gene giving information about the position in the genome (chr)
		 * For each perform a lookup in the fidx file (saved in "fasta_annotation" to gather the position more quickly.
		 *
		 * */
		try {


			RandomAccessFile raf = new RandomAccessFile(CommandLine_Parser.inputFile_fasta, "r");

			for (Map.Entry<String, Gen> gen : target_genes.entrySet()) {

/*
				System.out.println((gen.getValue().getStart()));
				System.out.println((gen.getValue().getEnd()));
*/

				gen.getValue().get_Transcripts().forEach((transcript_id, t) -> {
					StringBuilder stringBuilder = new StringBuilder();
/*
Loop ueber alle Exons eines Transkripts
*/
					t.getExons().forEach(exon -> {

						try {

//							raf.seek((fasta_annotation.get(gen.getValue().getchr())[1] + (exon.getStart() / fasta_annotation.get(gen.getValue().getchr())[2])) + exon.getStart());


							long lines_in_entry = (exon.getStart() - 1) / fasta_annotation.get(gen.getValue().getchr())[2];
							long last_line_length = exon.getStart() - (lines_in_entry * fasta_annotation.get(gen.getValue().getchr())[2]);
							long offset = lines_in_entry * fasta_annotation.get(gen.getValue().getchr())[3] + last_line_length;

							raf.seek(fasta_annotation.get(gen.getValue().getchr())[1] + offset - 1);


							long lines_in_1 = (int) Math.floor(exon.getEnd() / fasta_annotation.get(gen.getValue().getchr())[2]);
							long lines_in_2 = (int) Math.floor(exon.getStart() / fasta_annotation.get(gen.getValue().getchr())[2]);

							long lines_for_real = lines_in_1 - lines_in_2;

							long last_line = (exon.getEnd() - exon.getStart()) - (lines_for_real * fasta_annotation.get(gen.getValue().getchr())[2]);
							long off = lines_for_real * fasta_annotation.get(gen.getValue().getchr())[3] + last_line;


						/*	System.out.println(t.getTrans_id());
							System.out.println(t.getStart());
							System.out.println(t.getEnd());
							System.out.println("Exon:" + exon);
							System.out.println(fasta_annotation.get(gen.getValue().getchr())[1] + offset - 1);
							System.out.println("lines in " + lines_for_real);
							System.out.println("lastL:ine " + last_line);
							System.out.println("off " + off);
						*/
							int length = Math.toIntExact(off) + 1;
//							int length = Math.toIntExact(exon.getEnd()-exon.getStart()+1);


							byte[] bytey = new byte[length];


							raf.readFully(bytey, 0, length);
							String temp = new String(bytey, StandardCharsets.UTF_8).replaceAll("\n", "");
							stringBuilder.append(temp);

//TODO: Optimierung: Alles in einer for loop loesen.(generate reads and mutation)
						} catch (IOException e) {
							throw new TesException("Error Seeking correnct Line in FastaFile", e);
						}
					});

					t.add_Sequence(stringBuilder.toString());

/*
					System.out.println();
					System.out.println(stringBuilder.toString());
					System.out.println();
					System.out.println(t.get_Sequence());
*/

				});
			}

		} catch (IOException e) {
			throw new TesException("Failed to read Fasta File @RandomFileAccess", e);
		}

	}

	private void calculate_core() {

		/*
		 * TODO: 1. Berechnen der FL (Fragment length) mit Gaussian java utils
		 *
		 * Berechnen der FL (Fragment Length) mit Hilfe der gegebenen Parameter  (Standart derivate) und frlength.
		 * waehle eine zufaellige Position auf der Sequenz des Transkripts , die von 0 bis Transkript.length-FL liegt.
		 * Durch Readlength Parameter zwei Sequenzen generieren, jede so lang wie die Readlength. Diese koennen auch ueberlappen, muessen aber nicht.
		 * Die erste wird ab der Position des zufaelligen Pointer gelesen, die zweite vom Ende des Transkripts, dabei wird das reverse komplement gelesen, also fuer jede Base muss die komplementaere Base generiert werden.
		 * Anschliessend muessen fuer beide Sequenzen Mutationen berechnet werden.
		 * Die Positionen der Sequencen muessen anschliessend gespeichert werden, siehe hierzu ausgabe Files .
		 * Die Positionen muessen einmal in Genomische Daten umgerechnet werden. (da vorher nur Position relativ auf dem jeweiligen Chromosom angegeben wird).
		 * mapping.info enthaelt:
		 * ReadID   Chr Gene    Transkript  Transkriptpos.  t_fw_regvec fw_regvev   t_rw_regvec rw-regveg fw_mutations    rw_mutations
		 *readid	chr	gene	transcript	t_fw_regvec	t_rw_regvec	fw_regvec	rw_regvec	fw_mut	rw_mut
		 *
		 * Zu beachten: Sollte die FL kleiner als die Readlength sein, erneut berechnen.
		 * Simuliere Mutationen: Auslagern in extra Methode
		 *
		 * */
		int read_id = 0;
		Random r = new Random();

		for (Map.Entry<String, Gen> gen : target_genes.entrySet()) {

			for (Transcript t : gen.getValue().get_Transcripts().values()) {

				for (int i = readcounts.get(gen.getKey()).get(t.getTrans_id()) - 1; i >= 0; i--) {


					rand = new NormalDistribution(frlength, standardDeviation);


//					int fragmen_length = (int) Math.max(readlength, (r.nextGaussian() * standardDeviation + frlength));
					int fragmen_length;

					do {
						fragmen_length = (int) rand.sample();
					} while (fragmen_length <= readlength || fragmen_length + 1 >= t.get_Sequence().length());


					int random_pos = r.nextInt(t.get_Sequence().length() - fragmen_length);
					String sequence = t.get_Sequence(random_pos, fragmen_length);


					String fw = sequence.substring(0, readlength);
					String rw = reverse_komplement_calc(sequence.substring(sequence.length() - readlength, sequence.length()));

					String[] mutated_seq_fw = mutation_generator(read_id, fw);
					String[] mutated_seq_rw = mutation_generator(read_id, rw);

					fw_mutations_pointer.put(read_id, mutated_seq_fw[1]);
					rw_mutations_pointer.put(read_id, mutated_seq_rw[1]);

					fw_reads.put(read_id, mutated_seq_fw[0]);
					rw_reads.put(read_id, mutated_seq_rw[0]);

//TODO: Ausgeben der Coordinaten, lokal und chromosomal


					/*Auf Transkript*/
					int[] t_fw_regveg = {random_pos, random_pos + readlength};
					int[] t_rw_regveg = {random_pos + fragmen_length - readlength, random_pos + fragmen_length};

					/*Auf Gene/Chromosom*/
					String fw_regvec = t.get_Chromosomal_location(random_pos, random_pos + readlength);
					String rw_regvec = t.get_Chromosomal_location(random_pos + fragmen_length - readlength, random_pos + fragmen_length);


					System.out.println("Readlength: " + readlength);
					System.out.println("Transcript_id: " + t.getTrans_id());
					System.out.println("Sequence Length: " + t.get_Sequence().length());
//					System.out.println("Transcript Length: " + t.get_length());
					System.out.println("RandomPos: " + random_pos);
					System.out.println("fragment_length: " + fragmen_length);
//					System.out.println("Sequence: \n" + sequence);
//					System.out.println("forward read:\n" + fw);
//					System.out.println("forward read: Mutated \n" + mutated_seq_fw[0]);
//					System.out.println("backward read:\n" + rw);
//					System.out.println("backward read: Mutated \n" + mutated_seq_rw[0]);
					System.out.println();
					System.out.println("read_id " + read_id);
					System.out.println("fw_mutations_pointer: " + fw_mutations_pointer.get(read_id));
					System.out.println("rw_mutations_pointer: " + rw_mutations_pointer.get(read_id));
					System.out.println("fw_regvec: " + fw_regvec);
					System.out.println("rw_regvec: " + rw_regvec);


					System.out.println();

					read_id++;


				}
			}
		}
	}

	private String[] mutation_generator(int event_id, String sequence_neutral) {

		/*
		@return: [0]= mutated String, [1] = | seperated List of numbers, where mutations happened.
		Mutates String and returns String with mutations, adds readid_hashmap values
		* */

		StringBuilder mutationsocc = new StringBuilder();


		Random darwin = new Random();
		char[] beagle = sequence_neutral.toCharArray();


		int i = sequence_neutral.length() - 1;
		while (i >= 0) {
			if (darwin.nextFloat() <= mutationrate / 100) {

				int pointmutation = darwin.nextInt(sequence_neutral.length());
				int randomNum = ThreadLocalRandom.current().nextInt(1, 4);

				switch (sequence_neutral.charAt(i)) {
					case 'A':
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
								mutationsocc.append(i + ",");
								break;
							case 2:
								beagle[i] = 'G';
								mutationsocc.append(i + ",");
								break;
							case 3:
								beagle[i] = 'T';
								mutationsocc.append(i + ",");
								break;
						}
						break;
					case 'C':
						switch (randomNum) {
							case 1:
								beagle[i] = 'A';
								mutationsocc.append(i + ",");
								break;
							case 2:
								beagle[i] = 'G';
								mutationsocc.append(i + ",");
								break;
							case 3:
								beagle[i] = 'T';
								mutationsocc.append(i + ",");
								break;
						}
						break;
					case 'T':
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
								mutationsocc.append(i + ",");
								break;
							case 2:
								beagle[i] = 'G';
								mutationsocc.append(i + ",");
								break;
							case 3:
								beagle[i] = 'A';
								mutationsocc.append(i + ",");
								break;
						}
						break;
					case 'G':
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
								mutationsocc.append(i + ",");
								break;
							case 2:
								beagle[i] = 'A';
								mutationsocc.append(i + ",");
								break;
							case 3:
								beagle[i] = 'T';
								mutationsocc.append(i + ",");
								break;
						}
						break;
				}
			}
			i--;
		}
		if (mutationsocc.length() != 0)
			if (mutationsocc.charAt(mutationsocc.length() - 1) == ',') {
				mutationsocc.deleteCharAt(mutationsocc.length() - 1);
			}

		return new String[]{new String(beagle), mutationsocc.toString()};
	}

	private String reverse_komplement_calc(String sequence) {
		StringBuilder komp = new StringBuilder();
		for (int i = sequence.length() - 1; i >= 0; i--) {

			switch (sequence.charAt(i)) {

				case 'A':
					komp.append("T");
					break;
				case 'C':
					komp.append("G");
					break;
				case 'T':
					komp.append("A");
					break;
				case 'G':
					komp.append("C");
					break;
			}
		}
		return komp.toString();
	}

	private void printHeaders() {
//TODO: aufspalten in jeweiliges File
		/*
		 * To read.mappinginfo:
		 * */
		System.out.println("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regcev\tfw_regvec\trw_regveg\tfw_mut\trw_mut");

		/*
		 * to rw und fw Fasta: extra ausgabe. siehe print()
		 * */

	}


	public void print() {
//		TODO Set Filestreams and print to files

	}
}

