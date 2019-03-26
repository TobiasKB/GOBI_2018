import readsimulator.ReadSimulator;
import utils.FileHandlerUtil;

/**
 * @author TKB
 */
public class Runner {

	public static void main(String[] args) {

		FileHandlerUtil.parseParameters(args);
		/*
		try {

			PrintStream fileout = new PrintStream(FileHandlerUtil.outputFile);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		System.setOut(fileout);*/
//		Thread exonSkipping = new Thread(new Exon_Skipping_v2());
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
