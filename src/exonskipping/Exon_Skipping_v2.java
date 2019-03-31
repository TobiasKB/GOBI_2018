package exonskipping;

import utils.CommandLine_Parser;
import utils.GTF_FileParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class Exon_Skipping_v2 implements Runnable {


    public TreeMap<String, Gen> genes;
    public TreeMap<String, String> skippingEvents;


    public Exon_Skipping_v2() {
        genes = new TreeMap<String, Gen>();
        skippingEvents = new TreeMap<String, String>();
    }


    @Override
    public void run() {

        print_Output();

        build_Datastructure();

        collect_Skipping_Events();

        parse_Output();

    }


    public void build_Datastructure() {

        long speed_measure_start = System.currentTimeMillis();
        try {
            BufferedReader br = new BufferedReader(new FileReader(CommandLine_Parser.inputFile_GTF));
            String line = br.readLine();
            String[] cleanLine;

            while (line != null) {

                cleanLine = GTF_FileParser.parseLine(line, null, null);

                if (cleanLine != null) {
/**                    genes.put("sd",new Gen());
 cleanLine[0]= chr;
 cleanLine[1]= start;
 cleanLine[2]= end;
 cleanLine[3]= strand;
 cleanLine[4]= gene_name;
 cleanLine[5]= gene_id;
 cleanLine[6]= protein_id;
 cleanLine[7]= transcript_id;
 */
                    if (!(genes.containsKey(cleanLine[5]))) {
                        genes.put(cleanLine[5], new Gen(cleanLine[5], cleanLine[4], cleanLine[0], Integer.parseInt(cleanLine[1]), Integer.parseInt(cleanLine[2]), cleanLine[3], cleanLine[6]));
                    }

                    if (!genes.get(cleanLine[5]).get_Transcripts().containsKey(cleanLine[7]))
                        genes.get(cleanLine[5]).add_Transcript(new Transcript(cleanLine[7], cleanLine[3], cleanLine[6], Integer.parseInt(cleanLine[1]), Integer.parseInt(cleanLine[2]), cleanLine[5]));


                    genes.get(cleanLine[5]).add_Region(cleanLine[7], Integer.parseInt(cleanLine[1]), Integer.parseInt(cleanLine[2]), cleanLine[8]);
                    genes.get(cleanLine[5]).add_nprots();
                }
                line = br.readLine();
            }
        } catch (Exception e) {
            throw new RuntimeException("GTF File could not be read", e);
        }
        long execution_time = System.currentTimeMillis() - speed_measure_start;
//        System.out.printf("Took %.5f seconds", execution_time / 1e3);
    }

    public void collect_Skipping_Events() {

        genes.forEach((key, value) -> value.get_Transcripts().forEach((tkey, tentry) -> tentry.calculate_Introns()));
        genes.forEach((key, value) -> value.calculate_Skipping_Events());

    }

    public void parse_Output() {
        /*
         *       genID,
         *       gene_name,
         *       chr,
         *       strand,
         *       nprots,
         *       ntrans,
         *       SV,
         *       WT 1,
         *       WT 2,
         *       WT_prots,
         *       SV_prots,
         *       min_skipped_exon,
         *       max_skipped_exon,
         *       min_skipped_bases,
         *       max_skipped_bases
         * */


        StringBuilder stringBuilder = new StringBuilder();
        String output_line = "";
        List<String[]> out_lines = new ArrayList<>();
        for (Map.Entry<String, Gen> gene : genes.entrySet()) {

            String gene_id = gene.getKey();
            String gene_name = gene.getValue().getName();
            String chr = gene.getValue().getchr();
            String strand = gene.getValue().getstrand();

            String nprots = Integer.toString(gene.getValue().getnprot_size());
            String ntrans = Integer.toString(gene.getValue().count_trans);

            for (Map.Entry<String, List<String>> sv_to_trans : gene.getValue().sv_to_transcript_map.entrySet()) {
                String SV = sv_to_trans.getKey();
                String WT = gene.getValue().sv_to_wt_map.get(sv_to_trans.getKey()).toString().substring(1, gene.getValue().sv_to_wt_map.get(sv_to_trans.getKey()).toString().length() - 1).replaceAll("\\s*,\\s* ", "|");
                String WT_prots = "";
                String SV_prots = "";
                for (String svprot : gene.getValue().getSV_prots().get(sv_to_trans.getKey())) {
                    SV_prots += svprot;
                    SV_prots += "|";
                }
                SV_prots = SV_prots.substring(0, SV_prots.length() - 1);
                for (String stransid : gene.getValue().wt_prots.get(sv_to_trans.getKey())) {
                    WT_prots += stransid;
                    WT_prots += "|";
                }
                WT_prots = WT_prots.substring(0, WT_prots.length() - 1);
                String min_skipped_exon = Integer.toString(gene.getValue().sv_to_skippedExons.get(sv_to_trans.getKey())[1]);
                String max_skipped_exon = Integer.toString(gene.getValue().sv_to_skippedExons.get(sv_to_trans.getKey())[0]);
                String max_skipped_bases = Integer.toString(gene.getValue().sv_to_skippedBases.get(sv_to_trans.getKey())[0]);
                String min_skipped_bases = Integer.toString(gene.getValue().sv_to_skippedBases.get(sv_to_trans.getKey())[1]);

                System.out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", gene_id, gene_name, chr, strand, nprots, ntrans, SV, WT, WT_prots, SV_prots, min_skipped_exon, max_skipped_exon, min_skipped_bases, max_skipped_bases);
            }
        }
    }

    private void print_Output() {


        String header = "id\t" + "symbol\t" + "chr\t" + "strand\t" + "nprots\t" + "ntrans\t" + "SV\t" + "WT\t" + "WT_prots\t" + "SV_prots\t"
                + "min_skipped_exon\t" + "max_skipped_exon\t" + "min_skipped_bases\t" + "max_skipped_bases\t";
        System.out.println(header);
    }

    private void outPutParser() {
//genID, gene_name, chr, strand, nprots, ntrans, SV, WT, WT_prots, SV_prots, min_skipped_exon, max_skipped_exon, min_skipped_bases, max_skipped_bases
    }
}