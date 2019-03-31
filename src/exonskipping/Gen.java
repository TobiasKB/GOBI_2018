package exonskipping;

import java.util.*;

import static utils.HashUtil.mergeSet;

public class Gen {

    public int count_trans;
    public TreeMap<String, List<Intron>> sv_to_wt_map;
    public TreeMap<String, List<String>> sv_to_transcript_map;
    public HashMap<String, HashSet<String>> wt_prots;
    public HashMap<String, HashSet<String>> sv_prots;
    public HashMap<String, int[]> sv_to_skippedExons;
    public HashMap<String, int[]> sv_to_skippedBases;
    private String gene_id;
    private String chr;
    private String gene_name;
    private HashSet<String> proteins;
    private TreeMap<String, Transcript> tmap;
    private int ntrans;
    private int nprots;
    private int start;
    private int end;
    private String strand;

    public Gen(String gene_id, String gene_name, String chr, int start, int end, String strand, String protein_id) {
        tmap = new TreeMap<String, Transcript>();
        sv_to_transcript_map = new TreeMap<>();
        sv_to_wt_map = new TreeMap<>();
        proteins = new HashSet<String>();
        this.gene_name = gene_name;
        this.sv_prots = new HashMap<>();
        this.gene_id = gene_id;
        this.start = start;
        count_trans = 0;
        this.chr = chr;
        this.end = end;
        this.strand = strand;
        this.proteins.add(protein_id);
        this.wt_prots = new HashMap<>();
        this.sv_to_skippedExons = new HashMap<>();
        this.sv_to_skippedBases = new HashMap<>();

    }

    public void add_nprots() {
        this.nprots++;
    }

    public String getchr() {

        return this.chr;
    }

    public String getName() {
        return this.gene_name;
    }

    public String getID() {
        return this.gene_id;
    }

    public void add_Transcript(Transcript trans) {
        if (!(tmap.containsKey(trans.getTrans_id()))) {
            this.tmap.put(trans.getTrans_id(), trans);
            count_trans++;
        }
    }

    public void add_Region(String trans_id, int start, int end, String identifier) {
        tmap.get(trans_id).add_Region(start, end, identifier);
        update_end(end);
        update_start(start);
    }

    public void update_start(int new_start) {
        if (start > new_start)
            this.start = new_start;
    }

    public void update_end(int new_end) {
        if (end < new_end)
            this.end = new_end;
    }

    public TreeMap<String, Transcript> get_Transcripts() {
        return this.tmap;
    }

	public Transcript get_Transcript(String transcript_id) {
		return tmap.get(transcript_id);
	}

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public HashSet<String> get_proteins() {
        return this.proteins;
    }


    public void calculate_Skipping_Events() {

        tmap.forEach((key1, value1) -> tmap.forEach((key2, value2) -> value1.get_Interactions(value2)));

        tmap.forEach((transid, transcript) -> transcript.sv_to_wt_index.forEach((sv, wt_list) -> {
            if (!transcript.sv_to_wt_index.get(sv).isEmpty()) {
                if (!sv_to_skippedExons.containsKey(sv)) {
                    this.sv_to_skippedExons.put(sv, new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE});
                }

                this.sv_to_skippedExons.get(sv)[0] = Math.max(sv_to_skippedExons.get(sv)[0], transcript.sv_to_skippedExons.get(sv)[0]);

                this.sv_to_skippedExons.get(sv)[1] = Math.min(sv_to_skippedExons.get(sv)[1], transcript.sv_to_skippedExons.get(sv)[1]);

                if (!sv_to_skippedBases.containsKey(sv)) {
                    this.sv_to_skippedBases.put(sv, new int[]{Integer.MIN_VALUE, Integer.MAX_VALUE});
                }
                this.sv_to_skippedBases.get(sv)[0] = Math.max(sv_to_skippedBases.get(sv)[0], transcript.sv_to_skippedBases.get(sv)[0]);
                this.sv_to_skippedBases.get(sv)[1] = Math.min(sv_to_skippedBases.get(sv)[1], transcript.sv_to_skippedBases.get(sv)[1]);


                if (!wt_prots.containsKey(sv)) {
                    HashSet<String> wt_prots_temp = new HashSet<>();
                    wt_prots_temp = transcript.sv_to_wt_prots.get(sv);
                    wt_prots.put(sv, wt_prots_temp);
                } else {
                    wt_prots.put(sv, mergeSet(wt_prots.get(sv), transcript.sv_to_wt_prots.get(sv)));
                }

                if (!sv_to_wt_map.containsKey(sv)) {
                    sv_to_wt_map.put(sv, wt_list);
                }
                if (!sv_to_transcript_map.containsKey(sv)) {
                    List<String> sv_transcripts = new ArrayList<>();
                    sv_transcripts.add(transid);
                    sv_to_transcript_map.put(sv, sv_transcripts);

                } else {
                    sv_to_transcript_map.get(sv).add(transid);
                }
            }
        }));

        sv_to_transcript_map.forEach((sv, transcript_id) -> {
            transcript_id.forEach(x -> {
                if (!sv_prots.containsKey(sv)) {

                    HashSet<String> sv_prots_temp = new HashSet<>();
                    sv_prots_temp = get_Transcripts().get(x).get_proteins();

                    sv_prots.put(sv, sv_prots_temp);
                } else {
                    sv_prots.put(sv, mergeSet(sv_prots.get(sv), get_Transcripts().get(x).get_proteins()));
                }

            });
        });
    }

    public String getstrand() {
        return this.strand;
    }

    public int getnprot_size() {
        return this.nprots;

    }

    public String getSV_prots_toString() {
//        sv_prots.forEach((sv, prots)-> {
//            String svprots="";
//            prots.forEach(x->{svprots+= x;});
//            svprots+=prots.toString();
//        });
        return null;
    }


    public HashMap<String, HashSet<String>> getSV_prots() {
        return this.sv_prots;
    }


}