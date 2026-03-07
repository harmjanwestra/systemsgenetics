package transeqtlcrossmap;

import org.apache.commons.cli.*;

import java.io.IOException;

public class Main {


    public static void main(String[] args) {

        Options options = new Options();

        options.addRequiredOption("m", "mode", true, "Mode: [makeseq|quantify]");
        options.addRequiredOption("e", "eqtlfile", true, "eQTL file");
        options.addOption("g", "genomefile", true, "Genome fasta file");
        options.addOption("i", "indir", true, "Input directory");
        options.addRequiredOption("o", "out", true, "Output directory or file");
        options.addRequiredOption("a", "gtf", true, "GTF file");
        options.addOption("w", "refwindowsize", true, "Reference window size [default 5000000]");
        options.addOption("r", "readlen", true, "Read length [default 35]");
        options.addOption("s", "readshift", true, "Shift read window with [default 2]");
        options.addOption("t", "jobtemplate", true, "Job template file");
        options.addOption("c", "usecisgeneasrefwindow", false, "Use cis-gene as reference window [default: false]");

        TransEQTLCrossMapBulk mt = new TransEQTLCrossMapBulk();

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);


            String mode = cmd.getOptionValue("mode");

            System.out.println(mode);
            switch (mode) {
                case "makeseq":
                    int readlen = 35;
                    int readshift = 2;
                    int ciswindowsize = 5000000;
                    boolean exportindividualgenes = true;
                    boolean exportexonseqs = false;
                    boolean overwriteexistingsnps = true;

                    String genomefasta = null;
                    String gtfannotation = null;
                    String outdir = null;
                    String eqtlfile = null;
                    String jobtemplate = null;
                    boolean useCisGeneAsCisWindow = false;

                    boolean settingsOK = true;
                    if (cmd.hasOption("eqtlfile")) {
                        eqtlfile = cmd.getOptionValue("eqtlfile");
                    } else {
                        System.out.println("Required --eqtlfile");
                        settingsOK = false;
                    }
                    if (cmd.hasOption("genomefile")) {
                        genomefasta = cmd.getOptionValue("genomefile");
                    } else {
                        System.out.println("Required --genomefile");
                        settingsOK = false;
                    }
                    if (cmd.hasOption("out")) {
                        outdir = cmd.getOptionValue("out");
                    } else {
                        System.out.println("Required --out");
                        settingsOK = false;
                    }
                    if (cmd.hasOption("gtf")) {
                        gtfannotation = cmd.getOptionValue("gtf");
                    } else {
                        System.out.println("Required --gtf");
                        settingsOK = false;
                    }
                    if (cmd.hasOption("refwindowsize")) {
                        ciswindowsize = Integer.parseInt(cmd.getOptionValue("refwindowsize"));
                    } else {
                        System.out.println("Required --refwindowsize");
                        settingsOK = false;
                    }

                    if (cmd.hasOption("readlen")) {
                        readlen = Integer.parseInt(cmd.getOptionValue("readlen"));
                    }
                    if (cmd.hasOption("readshift")) {
                        readshift = Integer.parseInt(cmd.getOptionValue("readshift"));
                    }

                    if (cmd.hasOption("jobtemplate")) {
                        jobtemplate = cmd.getOptionValue("jobtemplate");
                    } else {
                        System.out.println("Required --jobtemplate");
                        settingsOK = false;
                    }
                    if (cmd.hasOption("usecisgeneasrefwindow")) {
                        useCisGeneAsCisWindow = true;
                    }

                    if (settingsOK) {
                        mt.createSequences(eqtlfile, genomefasta, gtfannotation, outdir, ciswindowsize, readlen, readshift,
                                exportindividualgenes, exportexonseqs, overwriteexistingsnps, jobtemplate, useCisGeneAsCisWindow);
                    }

                case "quantify":
                    outdir = null;
                    eqtlfile = null;
                    String indir = null;
                    gtfannotation = null;
                    settingsOK = true;
                    useCisGeneAsCisWindow = false;
                    ciswindowsize = 5000000;

                    if (cmd.hasOption("eqtlfile")) {
                        eqtlfile = cmd.getOptionValue("eqtlfile");
                    } else {
                        System.out.println("Required --eqtlfile");
                        settingsOK = false;
                    }

                    if (cmd.hasOption("out")) {
                        outdir = cmd.getOptionValue("out");
                    } else {
                        System.out.println("Required --out");
                        settingsOK = false;
                    }

                    if (cmd.hasOption("indir")) {
                        indir = cmd.getOptionValue("indir");
                    } else {
                        System.out.println("Required --indir");
                        settingsOK = false;
                    }

                    if (cmd.hasOption("gtf")) {
                        gtfannotation = cmd.getOptionValue("gtf");
                    } else {
                        System.out.println("Required --gtf");
                        settingsOK = false;
                    }

                    if (cmd.hasOption("usecisgeneasrefwindow")) {
                        useCisGeneAsCisWindow = true;
                    }

                    if (cmd.hasOption("refwindowsize")) {
                        ciswindowsize = Integer.parseInt(cmd.getOptionValue("refwindowsize"));
                    } else {
                        System.out.println("Required --refwindowsize");
                        settingsOK = false;
                    }

                    if (settingsOK) {
                        mt.quantifyIndividualEQTL(eqtlfile, indir, gtfannotation, outdir, useCisGeneAsCisWindow, ciswindowsize);
                    }
                default:
            }
        } catch (ParseException e) {
            System.out.println("Command line parse exception: " + e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.setWidth(160);
            formatter.printHelp("TransEQTLCrossMap.jar", options);
        } catch (IOException e) {
            System.out.println("Failed with IOException: " + e.getMessage());
            e.printStackTrace();
        } catch (Exception e) {
            System.out.println("Failed with Exception: " + e.getMessage());
            System.out.println("Probably some command line parameter was wrong. Try running java -jar MbQTL-SNAPSHOT-jar-with-dependencies.jar");
            e.printStackTrace();
        }
        System.out.println();
        System.out.println("Exiting...");
        System.out.println();


//        if (args.length < 1) {
//            System.out.println("Usage: makeseq || align || quantifyIndividualEQTL || quantifysnp || combine ");
//        } else if (args[0].equals("makeseq")) {
//
//            if (args.length < 10) {
//                System.out.println("Usage: makeseq eqtlfile genome ensembl.gtf.gz outdir windowsize[2000000] readlen[25] shift[2] exportindividualgenes exportexonseqs jobtemplate useCisGeneAsCisWindow");
//                System.out.println("Defaults: windowsize=2000000 readlen=25 shift=2");
//                System.exit(-1);
//            } else {
//                String eqtlfile = args[1]; // "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
//                String genomefasta = args[2]; // "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
//                String ensemblannotation = args[3]; // "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
//                String outdir = args[4]; // "D:\\Work\\eqtlgeneremap\\";
//                int windowsize = Integer.parseInt(args[5]); // 2000000;
//                int readlen = Integer.parseInt(args[6]); // 25;
//                int shift = Integer.parseInt(args[7]); // 2;
//                boolean exportindividualgenes = Boolean.parseBoolean(args[8]);
//                boolean exportexonseqs = Boolean.parseBoolean(args[9]);
//                boolean overwriteexistingsnps = true;
//                String jobtemplate = args[10];
//                boolean useCisGeneAsCisWindow = Boolean.parseBoolean(args[11]);
//                try {
//                    mt.createSequences(eqtlfile, genomefasta, ensemblannotation, outdir, windowsize, readlen, shift,
//                            exportindividualgenes, exportexonseqs, overwriteexistingsnps, jobtemplate, useCisGeneAsCisWindow);
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
//
////			String eqtlfile = "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
//////		String eqtlfile = "D:\\Work\\eqtlgeneremap\\test.txt";
////			String genomefasta = "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
////			String ensemblannotation = "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
////			String outdir = "D:\\Work\\eqtlgeneremap\\";
////			int windowsize = 2000000;
////			int readlen = 25;
////			int shift = 2;
//
//
//        } else if (args[0].equals("align")) {
//
//            if (args.length < 3) {
//                System.out.println("Usage: align indir eqtlfile");
//            } else {
//                String indir = args[1];
//                String eqtlfile = args[2];
//                try {
//                    mt.align(indir, eqtlfile);
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
//        } else if (args[0].equals("quantifyIndividualEQTL")) {
//            if (args.length < 5) {
//                System.out.println("Usage quantifyIndividualEQTL eqtlfile indir gtfannotation out");
//            } else {
//                String eqtlfile = args[1];
//                String indir = args[2];
//                String ensemblannotation = args[3];
//                String out = args[4];
//                try {
//                    mt.quantifyIndividualEQTL(eqtlfile, indir, ensemblannotation, out);
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
//        } else if (args[0].equals("quantifysnp")) {
//            if (args.length < 4) {
//                System.out.println("Usage quantifysnp eqtlfile indir out");
//            } else {
//                String alignment = args[1];
//                String genereadfile = args[2];
//                String out = args[3];
//                try {
//                    mt.quantifySNP(alignment, genereadfile, out);
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
//        } else if (args[0].equals("combine")) {
//
//            if (args.length < 4) {
//                System.out.println("Usage: combine alignout eqtl output");
//            } else {
//
//                String alignout = args[1];
//                String eqtl = args[2];
//                String output = args[3];
//                try {
//                    mt.combine(alignout, eqtl, output);
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
//            }
//
//        }


    }
}
