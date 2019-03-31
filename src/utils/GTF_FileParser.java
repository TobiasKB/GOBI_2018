package utils;

import java.util.StringTokenizer;
import java.util.TreeSet;

public final class GTF_FileParser extends Thread {
	/*
	 * File Parser for GTF Files
	 * @author TKB
	 **/
	private static String feature = "exon";
	private static TreeSet<Integer> tabs = new TreeSet<Integer>();
	private static String cdscheck, chr, strand, gene_id, protein_id, transcript_id, gene_name, start, end, exon_id = "";
	private static StringBuilder cleanliner = new StringBuilder();
	private static String[] cleanLine = new String[9];
	private static String filteroption_1;
	private static String filteroption_2;
	private static boolean filter1;
	private static boolean filter2;

	private GTF_FileParser() {
	}
//TODO: Rewrite method for "Exon" vs "CDS"
	/*
	 * Methode kopieren, eine fuer readSimulator, eine fuer ExonSkipping.
	 * Einmal mit  CDS , einmal mit Exon
	 * Exon: Hinzufuegen von Eigenschaften "exon_number" und "exon_id"
	 * ReadSimulator parseLine: Umstellen, Exon haben keine Protein_id --> im Transcript umstellen
	 *
	 * */

	public static String[] parseLine(String inline, String filteroption1, String filteroption2) {
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

		if (filteroption1 != null) {
			filteroption_1 = filteroption1;
			filter1 = true;
		}
		if (filteroption2 != null) {
			filteroption_2 = filteroption2;
			filter2 = true;
		}

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
				case "exon_id":
					exon_id = attr_splice.nextToken().replaceAll("\"", "");
//					System.out.println(exon_id);
					continue;
			}

		}
		/*
		 * Applies only for ReadSimulator TODO: REmove
		 * */
		if (protein_id.isEmpty()) {
			protein_id = "?";
		}
		if (exon_id.isEmpty()) {
			exon_id = "?";
		}
//        System.out.println();
		try {
			if (chr.isEmpty() | exon_id.isEmpty() | start.isEmpty() | end.isEmpty() | strand.isEmpty() | gene_id.isEmpty() | protein_id.isEmpty() | transcript_id.isEmpty() | gene_name.isEmpty()) {
				throw new IllegalStateException("Necessary Value for Datastructure is empty");
			}
		} catch (Exception e) {
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
		cleanLine[8] = exon_id;


		if (filter1) {
			if (!gene_id.equals(filteroption1)) {
				return null;
			}
		}

		if (filter2) {
			if (!transcript_id.equals(filteroption2)) {
				return null;
			}
		}
		return cleanLine;

//        return chr + "|" + start + "|" + end + "|" + strand + "|" + gene_id + "|" + protein_id + "|" + transcript_id + "|" + gene_name;

	}
}