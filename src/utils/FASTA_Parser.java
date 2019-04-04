package utils;

import java.util.HashMap;

public class FASTA_Parser {

	private static HashMap<String, long[]> fasta_annotation;

	//	bad practice
	private FASTA_Parser(HashMap fasta_annotation) {
		FASTA_Parser.fasta_annotation = fasta_annotation;

	}

	public static long[] newLine_koordinates(long start_position, long end_position, String chr) {

		long[] fasta_annotation_array = fasta_annotation.get(chr);


//		long lines_in_entry = ((start_position - 1) / fasta_annotation_array[2]);
//		long last_line_length = (start_position - (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]));
//		long offset =( ((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (start_position - (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2])));

		long start = fasta_annotation_array[1] + (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[3] + (start_position - (((start_position - 1) / fasta_annotation_array[2]) * fasta_annotation_array[2]))) - 1;


		long lines_in_1 = (int) Math.floor(end_position / fasta_annotation_array[2]);
		long lines_in_2 = (int) Math.floor(start_position / fasta_annotation_array[2]);

//		long lines_for_real = (lines_in_1 - lines_in_2);
//		long last_line = ((end_position - start_position) - ( (lines_in_1 - lines_in_2) * fasta_annotation_array[2]));
//		long off = ((lines_in_1 - lines_in_2) * fasta_annotation_array[3] + ((end_position - start_position) - ( (lines_in_1 - lines_in_2) * fasta_annotation_array[2])));

		long length = Math.toIntExact(((lines_in_1 - lines_in_2) * fasta_annotation_array[3] + ((end_position - start_position) - ((lines_in_1 - lines_in_2) * fasta_annotation_array[2])))) + 1;


		return new long[]{start, length};

	}

	public void setFasta_annotation() {

	}


}
