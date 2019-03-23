package utils;

/**
 * @author TKB
 */
public final class FileHandlerUtil {

    public static String inputFile;
    public static String outputFile;

    /**
     *
     */
    private FileHandlerUtil() {
        System.out.println("FileHandlerUtil cannot be invoked.");
    }


    public static void parseParameters(String... args) {
        check:
        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-gtf":
                    inputFile = args[i + 1];
//				System.out.println(inputFile);
                    i++;
                    break;
                case "-o":
                    outputFile = args[i + 1];
                    i++;
                    break;
                case "--998esf92mmow350sefmi048ssmlg4i":
                    System.out.println("You cracked it! ");
                    break;
                default:
                    System.out.println(
                            "Correct usage: -o \"output file path\" -g \"input file path\"");
                    break check;
            }
        }
    }


}
//-o "C:\Users\TKB\Downloads\Gobi_First_Assignment_Test\testoutput.txt" -gtf "C:\Users\TKB\Downloads\Gobi_First_Assignment_Test\Homo_sapiens.GRCh38.86.gtf"


//inputFile = "C:\\Users\\TKB\\Dropbox\\Studium\\Bioinformatik\\5.Semester\\Gobi\\ENSG00000131018_Unmod.txt";//
//inputFile="C:\\Users\\TKB\\Downloads\\Gobi_First_Assignment_Test\\Homo_sapiens.GRCh38.86.gtf";
