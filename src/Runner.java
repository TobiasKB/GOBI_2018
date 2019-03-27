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