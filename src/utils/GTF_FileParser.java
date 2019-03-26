package utils;

import java.util.StringTokenizer;
import java.util.TreeSet;

public final class GTF_FileParser extends Thread {
    /*
     * File Parser for GTF Files
     * @author TKB
     **/
    private static String feature = "CDS";
    private static TreeSet<Integer> tabs = new TreeSet<Integer>();
    private static String cdscheck, chr, strand, gene_id, protein_id, transcript_id, gene_name, start, end = "";
    private static StringBuilder cleanliner = new StringBuilder();
    private static String[] cleanLine = new String[8];

    private GTF_FileParser() {
    }


    public static String[] parseLine(String inline) {

        tabs.clear();
        cdscheck = "";
        chr = "";
        strand = "";
        gene_id = "";
        protein_id = "";
        transcript_id = "";
        gene_name = "";
        start = "";
        end = "";

        int ct = 0, hv1 = 0, c = 0;
        if (inline.charAt(0) == '#') {
            return null;
        }

        while (true) {
            hv1 = inline.indexOf('\t', ct);
            if (tabs.contains(hv1)) {
                break;
            }
            tabs.add(hv1);
            ct = hv1 + 1;
        }
        cdscheck = inline.substring(HashUtil.get_nthElement_Set(tabs, 2) + 1,
                HashUtil.get_nthElement_Set(tabs, 3));

        if (!(cdscheck.equals(feature))) {
            return null;
        }


        chr = inline.substring(0, HashUtil.get_nthElement_Set(tabs, 1));
        start = inline.substring(HashUtil.get_nthElement_Set(tabs, 3) + 1, HashUtil.get_nthElement_Set(tabs, 4));
        end = inline.substring(HashUtil.get_nthElement_Set(tabs, 4) + 1, HashUtil.get_nthElement_Set(tabs, 5));
        strand = inline.substring(HashUtil.get_nthElement_Set(tabs, 6) + 1, HashUtil.get_nthElement_Set(tabs, 7));

        String attributes = inline.substring(HashUtil.get_nthElement_Set(tabs, 8) + 1);
//        System.out.println(attributes);
        StringTokenizer attr_splice = new StringTokenizer(attributes, ";| ");
        while (attr_splice.hasMoreElements()) {

            switch (attr_splice.nextToken()) {
                case "gene_id":
                    gene_id = attr_splice.nextToken().replaceAll("\"", "");
//                    System.out.println(gene_id);
                    continue;
                case "protein_id":
                    protein_id = attr_splice.nextToken().replaceAll("\"", "");
//                    System.out.println(protein_id);
                    continue;
                case "transcript_id":
                    transcript_id = attr_splice.nextToken().replaceAll("\"", "");
//                    System.out.println(transcript_id);
                    continue;
                case "gene_name":
                    gene_name = attr_splice.nextToken().replaceAll("\"", "");
//                    System.out.println(gene_name);
                    continue;
            }

        }
//        System.out.println();
        try {
            if (chr.isEmpty() | start.isEmpty() | end.isEmpty() | strand.isEmpty() | gene_id.isEmpty() | protein_id.isEmpty() | transcript_id.isEmpty() | gene_name.isEmpty()) {
                throw new IllegalStateException("Necessary Value for Datastructure is empty");
            }
        } catch (IllegalStateException e) {
            throw new TesException("GTF_File_Parser ", e);
        }

        cleanLine[0] = chr;
        cleanLine[1] = start;
        cleanLine[2] = end;
        cleanLine[3] = strand;
        cleanLine[4] = gene_name;
        cleanLine[5] = gene_id;
        cleanLine[6] = protein_id;
        cleanLine[7] = transcript_id;
        return cleanLine;

//        return chr + "|" + start + "|" + end + "|" + strand + "|" + gene_id + "|" + protein_id + "|" + transcript_id + "|" + gene_name;

    }
}