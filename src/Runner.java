import exonskipping.Exon_Skipping_v2;
import utils.FileHandlerUtil;

import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author TKB
 */
public class Runner {

    /**
     * @param args
     * @throws InterruptedException
     */
    public static void main(String[] args) {
        FileHandlerUtil.parseParameters(args);
        PrintStream fileout = null;
        try {
            fileout = new PrintStream(FileHandlerUtil.outputFile);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        System.setOut(fileout);

        Thread t = new Thread(new Exon_Skipping_v2());
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
