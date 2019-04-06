package readsimulator;

import exonskipping.Gen;
import exonskipping.Transcript;
import org.apache.commons.math3.distribution.NormalDistribution;
import utils.CommandLine_Parser;
import utils.GTF_FileParser;
import utils.TesException;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
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
	private HashMap<Integer, String> read_mappinginfo;

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


	}

	private void readFiles() {
//		TODO: Stop everything until threads are done reading
		/*
		 * Already done by FileReader: get Paths and inputParameters
		 * Ein Thread liest das IndexFile , einer den Rest//
		 *
		 * */

	/*	ThreadPoolExecutor executor= (ThreadPoolExecutor) Executors.newFixedThreadPool(2);
		Task task_01 = new Task() {
			@Override
			protected Object call() throws Exception {
				return null;
			}
		};*/

		/*new Thread(new Runnable() {
			public void run() {*/
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
	/*		}
		}).start();

		new Thread(new Runnable() {
			public void run() {*/
		try {
//				Auslesen des readcount files
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
		/*	}
		}).start();
*/
		try {


			BufferedReader br = new BufferedReader(new FileReader(CommandLine_Parser.inputFile_GTF));
			String line = br.readLine();
			String[] cleanline;

			while (line != null) {

				cleanline = GTF_FileParser.parseLine(line, null, null);

				if (cleanline != null) {

//					for (String gene_id : readcounts.keySet()) {

//						for (String transcript_id : readcounts.get(gene_id).keySet()) {

//							if (cleanline[5].equals(gene_id) && cleanline[7].equals(transcript_id)) {

					if (!target_genes.containsKey(cleanline[5])) {

						target_genes.put(cleanline[5], new Gen(cleanline[5], cleanline[4], cleanline[0], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[3], cleanline[6]));

					}
					if (!target_genes.get(cleanline[5]).get_Transcripts().containsKey(cleanline[7]))
						target_genes.get(cleanline[5]).add_Transcript(new Transcript(cleanline[7], cleanline[3], cleanline[6], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[5], cleanline[0]));

					target_genes.get(cleanline[5]).add_Region(cleanline[7], Integer.parseInt(cleanline[1]), Integer.parseInt(cleanline[2]), cleanline[8]);
					target_genes.get(cleanline[5]).add_nprots();
//							}
//						}
//					}
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
			StringBuilder fileContent = new StringBuilder();


			FileWriter fw_fasta = new FileWriter(outputDirectory + "/fw.fastq", true);
			BufferedWriter fw_bf = new BufferedWriter(fw_fasta);
			PrintWriter fw_printer = new PrintWriter(fw_bf);

			FileWriter rw_fasta = new FileWriter(outputDirectory + "/rw.fastq", true);
			BufferedWriter rw_bf = new BufferedWriter(rw_fasta);
			PrintWriter rw_printer = new PrintWriter(rw_bf);

			FileWriter rmi = new FileWriter(outputDirectory + "/read.mappinginfo", true);
			BufferedWriter rmi_bf = new BufferedWriter(rmi);
			PrintWriter rmi_printer = new PrintWriter(rmi_bf);




			/*
			 *
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
			StringBuilder rudi = new StringBuilder();


			for (Gen gen : target_genes.values()) {

				if (!readcounts.containsKey(gen.getID())) {

					continue;

				}
				String fw, rw, sequence, fragment, fw_regvec, rw_regvec;
				int fragmen_length, random_pos;
				String[] mutated_seq_fw;
				String[] mutated_seq_rw;
				int[] t_fw_regveg = new int[2];
				int[] t_rw_regveg = new int[2];


				for (Transcript t : gen.get_Transcripts().values()) {

					if (!readcounts.get(gen.getID()).containsKey(t.getTrans_id())) {
						continue;

					}

					StringBuilder stringBuilder = new StringBuilder();
					String chr = gen.getchr();
					long[] fasta_annotation_array = fasta_annotation.get(chr);


					t.getExons().forEach(exon -> {
						long[] newLine_koordinates = new long[2];

						try {

//							raf.seek((fasta_annotation.get(gen.getValue().getchr())[1] + (exon.getStart() / fasta_annotation.get(gen.getValue().getchr())[2])) + exon.getStart());
							/*Calculates the correct location of the sequence and saves one sequences per transcript */

							newLine_koordinates = newLine_koordinates(exon.getStart(), exon.getEnd(), chr);

							raf.seek(newLine_koordinates[0]);


							byte[] bytey = new byte[(int) newLine_koordinates[1]];

							raf.readFully(bytey, 0, (int) newLine_koordinates[1]);
							String temp = new String(bytey, StandardCharsets.UTF_8).replaceAll("\n", "");
							stringBuilder.append(temp);

						} catch (IOException e) {
							throw new TesException("Error Seeking correnct Line in FastaFile", e);
						}
					});

					sequence = stringBuilder.toString();

					if (t.get_Strand().equals("-")) {
						sequence = reverse_komplement_calc(sequence);
					}

					t.add_Sequence(sequence);

					for (int i = readcounts.get(gen.getID()).get(t.getTrans_id()) - 1; i >= 0; i--) {

						rudi.delete(0, rudi.length());

						rand = new NormalDistribution(frlength, standardDeviation);

//					fragmen_length = (int) Math.max(readlength, (r.nextGaussian() * standardDeviation + frlength));

						do {
							fragmen_length = (int) rand.sample();
						} while (fragmen_length <= readlength || fragmen_length + 1 >= t.get_Sequence().length());

						random_pos = r.nextInt(t.get_Sequence().length() - fragmen_length);
						fragment = t.get_Sequence(random_pos, fragmen_length);


						fw = fragment.substring(0, readlength);
						rw = reverse_komplement_calc(fragment.substring(fragment.length() - readlength, fragment.length()));

						mutated_seq_fw = mutation_generator(read_id, fw);
						mutated_seq_rw = mutation_generator(read_id, rw);


						/*Auf Transkript*/
						t_fw_regveg[0] = random_pos;
						t_fw_regveg[1] = random_pos + readlength;

						t_rw_regveg[0] = random_pos + fragmen_length - readlength;
						t_rw_regveg[1] = random_pos + fragmen_length;

						/*Auf Gene/Chromosom*/

						fw_regvec = t.get_Chromosomal_location(random_pos, random_pos + readlength, fasta_annotation);
						rw_regvec = t.get_Chromosomal_location(random_pos + fragmen_length - readlength, random_pos + fragmen_length, fasta_annotation);


						if (fw_regvec.equals("") || rw_regvec.equals("")) {
							throw new TesException("Fehler beim regvev, gleich null");
						}
/*


						System.out.println("Transcript_id: " + t.getTrans_id()+"\t Transcript Strand: "+t.get_Strand());
						System.out.println("Transcript Exons: "+t.getExons());
						System.out.println("Sequenz of Transcript: \n" + t.get_Sequence());
						t.get_exon_to_local().forEach((exonid, startstop)->
										System.out.println(exonid+"\t"+startstop[0]+":"+startstop[1])
								);

						System.out.println("Readlength: " + readlength);
						System.out.println();
//						System.out.println("read_id " + read_id);
						System.out.println("Sequence Length: " + t.get_Sequence().length());
						System.out.println("Transcript Length: " + t.get_length());
						System.out.println("RandomPos: " + random_pos);
						System.out.println("fragment_length: " + fragmen_length);
						System.out.println("fragment: \n" + fragment);
//						System.out.println("forward read:\n" + fw);
						System.out.println("forward read: Mutated \n" + mutated_seq_fw[0]);
//						System.out.println("backward read:\n" + rw);
						System.out.println("backward read: Mutated \n" + mutated_seq_rw[0]);
						System.out.println();
						System.out.println("fw_mutations_pointer: " + mutated_seq_fw[1]);
						System.out.println("rw_mutations_pointer: " + mutated_seq_rw[1]);
						System.out.println("fw_regvec: " + fw_regvec);
						System.out.println("rw_regvec: " + rw_regvec);
						System.out.println();
*/

						//rudi causes too much Overhead. Kill rudi. You must.

						rudi.append(read_id + "\t");
						rudi.append(gen.getchr() + "\t");
						rudi.append(gen.getID() + "\t");
						rudi.append(t.getTrans_id() + "\t");
						rudi.append(t_fw_regveg[0] + "-" + t_fw_regveg[1] + "\t");
						rudi.append(t_rw_regveg[0] + "-" + t_rw_regveg[1] + "\t");
						rudi.append(fw_regvec + "\t");
						rudi.append(rw_regvec + "\t");
						rudi.append(mutated_seq_fw[1] + "\t");
						rudi.append(mutated_seq_rw[1] + "\t");
						rudi.append("\n");


			/*		System.out.println("Readlength: " + readlength);
					System.out.println("Transcript_id: " + t.getTrans_id());
					System.out.println("Sequence Length: " + t.get_Sequence().length());
					System.out.println("Transcript Length: " + t.get_length());
					System.out.println("RandomPos: " + random_pos);
					System.out.println("fragment_length: " + fragmen_length);
						System.out.println("Number of Exons in this transcript: "+t.getExons().size());
//					System.out.println("Sequence: \n" + t.get_Sequence());
//					System.out.println("forward read:\n" + fw);
//					System.out.println("forward read: Mutated \n" + mutated_seq_fw[0]);
//					System.out.println("backward read:\n" + rw);
//					System.out.println("backward read: Mutated \n" + mutated_seq_rw[0]);
//					System.out.println();
					System.out.println("read_id " + read_id);
					System.out.println("fw_mutations_pointer: " + mutated_seq_fw[1]);
					System.out.println("rw_mutations_pointer: " + mutated_seq_rw[1]);
					System.out.println("fw_regvec: " + fw_regvec);
					System.out.println("rw_regvec: " + rw_regvec);
					System.out.println();
*/


						fileContent.delete(0, fileContent.length());

						fileContent.append("@");
						fileContent.append(read_id + "\n");
						fileContent.append(mutated_seq_fw[0] + "\n");
						fileContent.append("+" + read_id + "\n");
						for (int z = 0; z < mutated_seq_fw[0].length(); z++) {
							fileContent.append("I");
						}
						fileContent.append("\n");


						fw_printer.print(fileContent);


						fileContent.delete(0, fileContent.length());


						fileContent.append("@");
						fileContent.append(read_id + "\n");
						fileContent.append(mutated_seq_rw[0] + "\n");
						fileContent.append("+" + read_id + "\n");
						for (int x = 0; x < mutated_seq_rw[0].length(); x++) {
							fileContent.append("I");
						}
						fileContent.append("\n");


						rw_printer.print(fileContent);

						fileContent.delete(0, fileContent.length());


						fileContent.append(rudi.toString());


						rmi_printer.print(fileContent);


						read_id++;


					}
				}

			}


			fw_printer.close();
			rw_printer.close();
			rmi_printer.close();

			fw_bf.close();
			rw_bf.close();
			rmi_bf.close();

			fw_fasta.close();
			rw_fasta.close();
			rmi.close();


		} catch (IOException e) {
			throw new TesException("Failed to read Fasta File @RandomFileAccess", e);
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

	private long[] newLine_koordinates(long start_position, long end_position, String chr) {

		long[] fasta_annotation_array = fasta_annotation.get(chr);

		long start = fasta_annotation_array[1] + (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (start_position - (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;

		long length = Math.toIntExact(((((int) Math.floor(end_position / fasta_annotation_array[2])) - ((int) Math.floor(start_position / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((end_position - start_position) - ((((int) Math.floor(end_position / fasta_annotation_array[2])) - ((int) Math.floor(start_position / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;

		return new long[]{start, length};

	}

	private void printHeaders() {
		try {

			RandomAccessFile print_stream_read_mappinfo = new RandomAccessFile(outputDirectory + "/read.mappinginfo", "rw");
			FileChannel channel_3 = print_stream_read_mappinfo.getChannel();
			channel_3.position(channel_3.size());
			StringBuilder fileContent = new StringBuilder();

			fileContent.append("readid\tchr_id\tgene_id\ttranscript_id\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut\n");

			byte[] strBytes = fileContent.toString().getBytes();
			ByteBuffer buffy = ByteBuffer.allocate(strBytes.length);
			buffy.put(strBytes);
			buffy.flip();
			channel_3.write(buffy);
			print_stream_read_mappinfo.close();
			channel_3.close();


		} catch (IOException e) {
			throw new TesException("Could not Read inputFile_fidx", e);
		}
	}
}

