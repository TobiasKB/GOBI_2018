package exonskipping;

import utils.TesException;

import java.util.*;

import static utils.HashUtil.mergeSet;

public class Transcript {


	public TreeMap<String, List<Intron>> sv_to_wt_index;
	public HashMap<String, int[]> sv_to_skippedBases;
	public HashMap<String, int[]> sv_to_skippedExons;
	public HashMap<String, HashSet<String>> sv_to_wt_prots;
	public RegionVector transvec;
	private String trans_id;
	private String strand;
	private String source;
	private TreeSet<Exon> exons;
	private TreeSet<Intron> introns;
	private HashSet<String> proteins;
	private int start;
	private int end;
	private String sequence;
	private int length_exons;
	private int real_length;
	private String chr;

	/*
	 * @exon_to_localMap: Stores the Exon_id and maps it to the lokal coordinates of the exon (relative to length_exons)
	 * */
	private LinkedHashMap<String, int[]> exon_to_lokalMap;


	public Transcript(String trans_id, String strand, String proteinID, int start, int end, String source, String chr) {
		this.trans_id = trans_id;
//		System.out.println("erschaffe Transkript mit id: "+trans_id);
		this.strand = strand;
		this.chr = chr;
		/*if(strand.equals("-")){
			exons = new TreeSet<>(Collections.reverseOrder());
		}else{
		}*/
		exons = new TreeSet<>();
		introns = new TreeSet<Intron>();
		this.start = start;
		this.end = end;
		this.proteins = new HashSet<String>();
		this.sv_to_wt_index = new TreeMap<String, List<Intron>>();
		this.sv_to_skippedBases = new HashMap<>();
		this.sv_to_skippedExons = new HashMap<>();
		this.sv_to_wt_prots = new HashMap<>();
		this.proteins.add(proteinID);
		this.source = source;
		this.length_exons = 1;
		this.real_length = 0;
		this.exon_to_lokalMap = new LinkedHashMap<String, int[]>();

	}


	public LinkedHashMap<String, int[]> get_exon_to_local() {
		return this.exon_to_lokalMap;
	}

	public void add_Sequence(String sequence) {
		this.sequence = sequence;
	}

	public String get_Sequence() {
		return this.sequence;
	}

	public String get_Sequence(int position, int length) {
		return sequence.substring(position, position + length);
	}

	public HashSet<String> get_proteins() {
		return this.proteins;
	}

	public void add_Region(int start, int stop, String identifier) {

		exons.add(new Exon(start, stop, identifier));
		if (start < this.start) {
			this.start = start;
		}
		if (stop > this.end) {
			this.end = stop;
		}
//		System.out.println(trans_id+" neues Exon: "+start+":"+stop);
//		exon_to_lokalMap.put(identifier, new int[]{this.length_exons, length_exons + (stop - start) + 1});
		/*System.out.println(this.length_exons);
		System.out.println(stop-start);


		*/
		exon_to_local_sort();

//		this.length_exons = +this.length_exons + stop - start + 1;
		this.real_length = this.end - this.start;

	}

	private void exon_to_local_sort() {

		this.length_exons = 0;
		for (Exon ex : exons) {
			exon_to_lokalMap.put(ex.get_ID(), new int[]{this.length_exons, length_exons + (ex.getEnd() - ex.getStart()) + 1});
			this.length_exons = this.length_exons + ex.getEnd() - ex.getStart() + 1;
		}

		int exTemplength = this.length_exons;
		if (this.strand.equals("-")) {
			for (Exon ex : exons) {

				exon_to_lokalMap.put(ex.get_ID(), new int[]{exTemplength - (ex.getEnd() - ex.getStart()) - 1, exTemplength});
				exTemplength = exTemplength - (ex.getEnd() - ex.getStart()) - 1;

			}
		}
	}


	public String get_Strand() {
		return this.strand;
	}

	public String get_Chromosomal_location(int local_start, int local_stop, HashMap<String, long[]> fasta_annotation) {
		System.out.println("Berechne Chromosomale Location");

		int[] koordinate = new int[2];
		StringBuilder koordinates = new StringBuilder();

		long[] fasta_annotation_array = fasta_annotation.get(this.chr);

		for (Exon ex : exons) {
			System.out.println("Transcript_id: " + trans_id + " Transcript Strand: " + strand);
			System.out.println(exons);
			System.out.println(ex.get_ID());
			System.out.println(ex);
			System.out.println("localstart:" + local_start);
			System.out.println("local stop:" + local_stop);
			System.out.println("Current Exon Start: "+exon_to_lokalMap.get(ex.get_ID())[0]);
			System.out.println("Current Exon Stop: "+exon_to_lokalMap.get(ex.get_ID())[1]);
			System.out.println("Exon Genomic Start: " + ex.getStart());
			System.out.println("Exon Genomic Stopp: " + ex.getEnd());
//			Ende des Exons bereits vor Start
			if (exon_to_lokalMap.get(ex.get_ID())[1] < local_start) {
				continue;
			}
//			Ende lokal bereits vor Exon Anfang
			if (exon_to_lokalMap.get(ex.get_ID())[0] > local_stop) {
//				return koordinates.toString();
				continue;
			}


			if (exon_to_lokalMap.get(ex.get_ID())[0] < local_start) {

				if (exon_to_lokalMap.get(ex.get_ID())[1] >= local_stop) {
					koordinate[0] = ex.getStart() + (local_start - exon_to_lokalMap.get(ex.get_ID())[0]);
					koordinate[1] = koordinate[0] + (local_stop - local_start);
					/*long start = fasta_annotation_array[1] + (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (koordinate[0] - (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;
					long length = Math.toIntExact(((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((koordinate[1] - koordinate[0]) - ((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;
					long end = start+length;*/
					koordinates.append(koordinate[0] + "-" + koordinate[1]);
//					koordinates.append(start + "-" + end);
//					System.out.println(koordinates.toString());
					return koordinates.toString();

				} else {
					koordinate[0] = ex.getEnd() - ((exon_to_lokalMap.get(ex.get_ID())[1]) - local_start);
					koordinate[1] = ex.getEnd();
					local_start = local_start + (koordinate[1] - koordinate[0]);//schiebt local_start auf start des naechsten Exons
					/*long start = fasta_annotation_array[1] + (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (koordinate[0] - (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;
					long length = Math.toIntExact(((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((koordinate[1] - koordinate[0]) - ((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;
					long end = start+length;
*/
					koordinates.append(koordinate[0] + "-" + koordinate[1] + "|");
//					koordinates.append(start + "-" + end + "|");
				}


			} else {

				if (exon_to_lokalMap.get(ex.get_ID())[1] > local_stop) {

					koordinate[0] = ex.getStart();
					koordinate[1] = koordinate[0] + (local_stop - local_start);
					/*long start = fasta_annotation_array[1] + (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (koordinate[0] - (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;
					long length = Math.toIntExact(((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((koordinate[1] - koordinate[0]) - ((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;
					long end = start+length;*/
					koordinates.append(koordinate[0] + "-" + koordinate[1]);
//					koordinates.append(start + "-" + end);
//					System.out.println(koordinates.toString());
					return koordinates.toString();

				} else if (exon_to_lokalMap.get(ex.get_ID())[1] == local_stop) {

					koordinate[0] = ex.getStart();
					koordinate[1] = ex.getEnd();/*
					long start = fasta_annotation_array[1] + (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (koordinate[0] - (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;
					long length = Math.toIntExact(((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((koordinate[1] - koordinate[0]) - ((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;
					long end = start+length;*/
					koordinates.append(koordinate[0] + "-" + koordinate[1]);
//					koordinates.append(start + "-" + end);
//					System.out.println(koordinates.toString());
					return koordinates.toString();


				} else if (exon_to_lokalMap.get(ex.get_ID())[1] < local_stop) {
					koordinate[0] = ex.getStart();
					koordinate[1] = koordinate[0] + (local_stop - local_start);
					local_start = local_start + (koordinate[1] - koordinate[0]);
					/*long start = fasta_annotation_array[1] + (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (koordinate[0] - (((koordinate[0] - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;
					long length = Math.toIntExact(((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[3] + ((koordinate[1] - koordinate[0]) - ((( (int) Math.floor(koordinate[1] / fasta_annotation_array[2])) - ((int) Math.floor(koordinate[0] / fasta_annotation_array[2]))) * fasta_annotation_array[2])))) + 1;
					long end = start+length;*/
					koordinates.append(koordinate[0] + "-" + koordinate[1] + "|");
//					koordinates.append(start + "-" + end + "|");
				}


			}


		}
		return "WRONG WRONG WRONG WRONG";
	}


	public void calculate_Introns() {
		Region temp = null;
		for (Exon ex : exons) {
			if (temp != null) {
				introns.add(new Intron((temp.end + 1), ex.start));
			}
			temp = ex;
		}

	}

	public String getTrans_id() {
		return trans_id;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public void get_Interactions(Transcript t2) {

		TreeMap<String, List<Intron>> sv_wt_map_temp = new TreeMap<>();

		ArrayList<Intron> intron_list_temp = new ArrayList<>();

		if (this.trans_id.equals(t2.trans_id)) {
			return;
		}
		if (this.start > t2.getEnd() || this.end < t2.getStart()) {
			return;
		}

		for (Intron intron : introns) {

			for (Intron intron_other : t2.get_Introns()) {

				if (intron.start == intron_other.start || intron.end == intron_other.end) {

					if (intron.equals(intron_other)) {
						continue;
					}
					if (intron.start < intron_other.start || intron.end > intron_other.end) {
						if (!sv_wt_map_temp.containsKey(intron.toString())) {
							intron_list_temp = new ArrayList<>();
							intron_list_temp.add(intron_other);
							sv_wt_map_temp.put(intron.toString(), intron_list_temp);
						}
						if (!intron_list_temp.contains(intron_other)) {
							intron_list_temp.add(intron_other);
							sv_wt_map_temp.get(intron.toString()).add(intron_other);

						}
					}
				}
			}
			if (sv_wt_map_temp.containsKey(intron.toString())) {

				if (sv_wt_map_temp.get(intron.toString()).size() > 1) {

					t2.get_Introns().forEach(x -> {

						if (intron.start < x.start && intron.end > x.end) {
							if (sv_wt_map_temp.containsKey(intron.toString())) {
								sv_wt_map_temp.get(intron.toString()).add(x);
							}
						}

					});


             /*       if (!sv_to_skippedBases.containsKey(intron.toString())) {
                        sv_to_skippedBases.put(intron.toString(), new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE});
                    }


                    int skipped_bases = calculate_skipps(sv_wt_map_temp.get(intron.toString()));

                    sv_to_skippedBases.get(intron.toString())[0] = Math.max(sv_to_skippedBases.get(intron.toString())[0], skipped_bases);
                    sv_to_skippedBases.get(intron.toString())[1] = Math.min(sv_to_skippedBases.get(intron.toString())[1], skipped_bases);

*/


					if (!sv_to_wt_prots.containsKey(intron.toString())) {
						HashSet<String> wt_prots_temp = new HashSet<>();
						wt_prots_temp = t2.proteins;
						sv_to_wt_prots.put(intron.toString(), wt_prots_temp);
					} else {
						sv_to_wt_prots.put(intron.toString(), mergeSet(t2.proteins, sv_to_wt_prots.get(intron.toString())));
					}

					if (!sv_to_wt_index.containsKey(intron.toString())) {
						sv_to_wt_index.put(intron.toString(), sv_wt_map_temp.get(intron.toString()));

					} else {
						sv_wt_map_temp.get(intron.toString()).removeAll(sv_to_wt_index.get(intron.toString()));
						sv_to_wt_index.get(intron.toString()).addAll(sv_wt_map_temp.get(intron.toString()));
					}
					try {
						for (Exon exon_other : t2.getExons()) {
							if (exon_other.start > intron.start && exon_other.end < intron.end) {
								int temp_ex_size = new RegionVector(this.exons, intron.start, intron.end).size;
								int temp_ex_o_size = new RegionVector(t2.getExons(), intron.start, intron.end).size;

								if (temp_ex_size <= temp_ex_o_size) {

									if (!sv_to_skippedExons.containsKey(intron.toString())) {
										sv_to_skippedExons.put(intron.toString(), new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE});
									}

									int max_Exons = Math.max(temp_ex_o_size - temp_ex_size, sv_to_skippedExons.get(intron.toString())[0]);
									int min_Exons = Math.min(temp_ex_o_size - temp_ex_size, sv_to_skippedExons.get(intron.toString())[1]);
									sv_to_skippedExons.get(intron.toString())[0] = Math.max(max_Exons, 1);
									sv_to_skippedExons.get(intron.toString())[1] = Math.max(min_Exons, 1);

									if (!sv_to_skippedBases.containsKey(intron.toString())) {
										sv_to_skippedBases.put(intron.toString(), new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE});
									}
									int skipped_bases = new RegionVector(t2.getExons(), intron.start, intron.end).length - new RegionVector(this.exons, intron.start, intron.end).length;
									sv_to_skippedBases.get(intron.toString())[0] = Math.max(sv_to_skippedBases.get(intron.toString())[0], skipped_bases);
									sv_to_skippedBases.get(intron.toString())[1] = Math.min(sv_to_skippedBases.get(intron.toString())[1], skipped_bases);
								}
							}
						}
					} catch (Exception e) {
						new TesException("Fehler beim Berechnen der Skipps", e);
					}
				}
			}
		}
		sv_wt_map_temp.clear();
	}

	public TreeSet<Exon> getExons() {
		return this.exons;
	}

	public void printEvents() {
		System.out.println("Transcript_id: " + trans_id + "\n" + sv_to_wt_index + "");
		System.out.println(sv_to_skippedBases.keySet());
		System.out.println();
	}

	public void print_Regions(String id) {
		System.out.println("Transkript:" + id + " \nIntrons:");
		for (Intron i : introns) {
			System.out.println(i);
		}
		System.out.println("Exons:");
		for (Exon e : exons) {
			System.out.println(e);
		}
		System.out.println();
	}

	public TreeSet<Intron> get_Introns() {
		return this.introns;
	}

	public int get_length() {
		return this.length_exons;
	}

	public int get_realLength() {
		return this.real_length;
	}


}