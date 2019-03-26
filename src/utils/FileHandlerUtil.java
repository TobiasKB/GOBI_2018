package utils;

/**
 * @author TKB
 */
public final class FileHandlerUtil {

    public static String inputFile_GTF;
    public static String inputFile_fasta;
    public static String inputFile_fidx;
    public static String inputFile_readcounts;
    public static String outputFile;
    public static String outputDirectory;
    public static int readlength;
    public static int frlength;
    public static int standardDeviation;
    public static double mutationrate;

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
                    inputFile_GTF = args[i + 1];
//				System.out.println(inputFile);
                    i++;
                    break;
                case "-o":
                    outputFile = args[i + 1];
                    i++;
                    break;
                case "-od":
                    outputDirectory = args[i + 1];
                    i++;
                    break;
                case "-fasta":
                    inputFile_fasta = args[i + 1];
                    i++;
                    break;
                case "-fidx":
                    inputFile_fidx = args[i + 1];
                    i++;
                    break;
                case "-readcounts":
                    inputFile_readcounts = args[i + 1];
                    i++;
                    break;
                case "-length":
                    readlength = Integer.parseInt(args[i + 1]);
                    i++;
                    break;

                case "-mutationrate":
                    mutationrate = Double.parseDouble(args[i + 1]);
                    i++;
                    break;
                case "-frlength":
                    frlength = Integer.parseInt(args[i + 1]);
                    i++;
                    break;
                case "-SD":
                    standardDeviation = Integer.parseInt(args[i + 1]);
                    i++;
                    break;
                case "--998esf92mmow350sefmi048ssmlg4i":
                    System.out.println("You cracked it! ");
                    break;
                default:
                    System.out.println(
//                            "Correct usage: -o \"output file path\" -g \"input file path\"");
                            "Correct usage: -o \"output file path\" -g \"input file path\" -od \"output directory path\" -gtf \"input file path to GTF File\" -fasta \"input file path to fasta\" " +
                                    "-mutationrate \"mutationrate\" -fidx \"input file path to fasta annotation file\" -readcounts \"input file path to readcounts\" -frlength \"fragment length\" -SD \"Standard Derivation\" -length \"Readlengths\"  ");
                    break check;
            }
        }
    }


}
//-o "C:\Users\TKB\Downloads\Gobi_First_Assignment_Test\testoutput.txt" -gtf "C:\Users\TKB\Downloads\Gobi_First_Assignment_Test\Homo_sapiens.GRCh38.86.gtf"


//inputFile = "C:\\Users\\TKB\\Dropbox\\Studium\\Bioinformatik\\5.Semester\\Gobi\\ENSG00000131018_Unmod.txt";//
//inputFile="C:\\Users\\TKB\\Downloads\\Gobi_First_Assignment_Test\\Homo_sapiens.GRCh38.86.gtf";
