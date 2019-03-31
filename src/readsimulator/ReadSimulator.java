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
	private HashMap<String, String> fw_reads;
	private HashMap<String, String> rw_reads;

	/*
	 * @fw_mutations:<read_id,  mutations>
	 * @rw_mutations:<read_id,  mutations>
	 * */
	private HashMap<String, int[]> fw_mutations;
	private HashMap<String, int[]> rw_mutations;


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

//TODO to shorten runtime, call the fragment Length Calculation right here and do not specificlly call it twice /
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

	private void calculateFragments() {

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

		Random r = new Random();
		for (Map.Entry<String, Gen> gen : target_genes.entrySet()) {
			gen.getValue().get_Transcripts().forEach((transcript_id, t) -> {
//TODO: Transcript length != sequence length!! --> get sequence length for now ?
				int fragmen_length = (int) Math.max(readlength, (r.nextGaussian() * standardDeviation + frlength));
				int random_pos = r.nextInt(t.get_Sequence().length() - fragmen_length);
				String sequence = t.get_Sequence(random_pos, fragmen_length);


				String fw = sequence.substring(0, readlength);
				String rw = reverse_komplement_calc(sequence.substring(sequence.length() - readlength, sequence.length()));
				/*System.out.println(t.getTrans_id());
				System.out.println("Sequence Length: "+t.get_Sequence().length());
				System.out.println("Transcript Length: "+t.get_length());
				System.out.println("RandomPos: "+ random_pos);
				System.out.println("fragment_length: "+ fragmen_length);
				System.out.println("Sequence: \n"+sequence);
				System.out.println("forward read:\n"+fw);
				System.out.println("backward read:\n"+rw);
				System.out.println();*/

//TODO: Entwerder direkt printen oder aber erst noch in HashMap abspeichern. Speichern in Hash!
//			TODO: Generate mutations
			});
		}
	}


	private HashMap<String, int[]> mutation_generator(String event_id, String sequence_neutral) {

		/*

		Return:
		1. Mutierte Sequenz>> in Hashmap abspeichern
		2.

		* HashMap mutation_list: <read_id, [fw_mutations_index],[rw_mutations_index]
		* problem:  Ein Read kann x viele mutationen enthalen --> List<int>
		* Was muss gespeichert werden ?
		* 1. wo
		* "Was" muss nicht gespeichert werden .
		*
		* <String, fw_list, rw_list >
		* */
		HashMap<String, int[]> mutation_list = new HashMap<>();


		Random darwin = new Random();
		char[] beagle = sequence_neutral.toCharArray();


		int i = sequence_neutral.length() - 1;
		while (i >= 0) {
			if (darwin.nextFloat() <= mutationrate) {

				int pointmutation = darwin.nextInt(sequence_neutral.length());
				int randomNum = ThreadLocalRandom.current().nextInt(1, 3);

				System.out.println("Random Number with ThreadLovalRandom: " + randomNum);

				switch (sequence_neutral.charAt(i)) {
					case 'A':
//						TODO: Erzeuge Random zwischen 1 und 3 => ueberpruefen!
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
							case 2:
								beagle[i] = 'G';
							case 3:
								beagle[i] = 'T';
						}
						break;
					case 'C':
						switch (randomNum) {
							case 1:
								beagle[i] = 'A';
							case 2:
								beagle[i] = 'G';
							case 3:
								beagle[i] = 'T';
						}
						break;
					case 'T':
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
							case 2:
								beagle[i] = 'G';
							case 3:
								beagle[i] = 'A';
						}
						break;
					case 'G':
						switch (randomNum) {
							case 1:
								beagle[i] = 'C';
							case 2:
								beagle[i] = 'A';
							case 3:
								beagle[i] = 'T';
						}
						break;
				}
			}
			i--;
		}
		return null;
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

//		TODO: Print Headers to outputfiles; generate outputfiles. Remember to concat, not overwrite them later on!


	}


}

