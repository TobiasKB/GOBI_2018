import readsimulator.ReadSimulator;
import utils.CommandLine_Parser;

import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author TKB
 */
public class Runner {

	public static void main(String[] args) {
		CommandLine_Parser.parseParameters(args);

		try {
			PrintStream fileout = new PrintStream(CommandLine_Parser.outputFile);
			System.setOut(fileout);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
//		Thread t = new Thread(new Exon_Skipping_v2());
		Thread t = new Thread(new ReadSimulator());
		t.start();
		synchronized (t) {
			try {
				t.wait();
			} catch (Exception e) {
				throw new RuntimeException("Startup of Program Failure", e);
			}
		}
	}

}
//"/home/tobias/Documents/Unterlagen Universität/GOBI_2018/Assignment_02/ReadSimulator/Test_Output/output"
// 	-length 75 -frlength 20 -SD 80 -mutationrate 1.0 -gtf /home/k/kraftblank/Documents/5.Semester/GOBI_2018_Sources/ReadSimulator/Homo_sapiens.GRCh37.75.gtf -fasta /home/k/kraftblank/Documents/5.Semester/GOBI_2018_Sources/ReadSimulator/Homo_sapiens.GRCh37.75.dna.toplevel.fa -fidx /home/k/kraftblank/Documents/5.Semester/GOBI_2018_Sources/ReadSimulator/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai -readcounts /home/k/kraftblank/Documents/5.Semester/GOBI_2018_Sources/ReadSimulator/simul.readcons -od /home/k/kraftblank/Documents/5.Semester/GOBI_2018_Sources/output.testfile
