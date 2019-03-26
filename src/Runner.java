import exonskipping.Exon_Skipping_v2;
import utils.CommandLine_Parser;

/**
 * @author TKB
 */
public class Runner {

	public static void main(String[] args) {
		CommandLine_Parser.parseParameters(args);

		Thread t = new Thread(new Exon_Skipping_v2());
//		Thread t = new Thread(new ReadSimulator());
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
