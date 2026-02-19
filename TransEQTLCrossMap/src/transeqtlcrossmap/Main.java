package transeqtlcrossmap;

import java.io.IOException;

public class Main {

    public static void main(String[] args){

        TransEQTLCrossMapBulk mt = new TransEQTLCrossMapBulk();
        if (args.length < 1) {
            System.out.println("Usage: makeseq || align || quantifyIndividualEQTL || quantifysnp || combine ");
        } else if (args[0].equals("makeseq")) {

            if (args.length < 10) {
                System.out.println("Usage: makeseq eqtlfile genome ensembl.gtf.gz outdir windowsize[2000000] readlen[25] shift[2] exportindividualgenes exportexonseqs jobtemplate");
                System.out.println("Defaults: windowsize=2000000 readlen=25 shift=2");
                System.exit(-1);
            } else {
                String eqtlfile = args[1]; // "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
                String genomefasta = args[2]; // "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
                String ensemblannotation = args[3]; // "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
                String outdir = args[4]; // "D:\\Work\\eqtlgeneremap\\";
                int windowsize = Integer.parseInt(args[5]); // 2000000;
                int readlen = Integer.parseInt(args[6]); // 25;
                int shift = Integer.parseInt(args[7]); // 2;
                boolean exportindividualgenes = Boolean.parseBoolean(args[8]);
                boolean exportexonseqs = Boolean.parseBoolean(args[9]);
                boolean overwriteexistingsnps = true;
                String jobtemplate = args[10];
                try {
                    mt.createSequences(eqtlfile, genomefasta, ensemblannotation, outdir, windowsize, readlen, shift,
                            exportindividualgenes, exportexonseqs, overwriteexistingsnps, jobtemplate);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

//			String eqtlfile = "D:\\Work\\eqtlgeneremap\\eQTLs-allsnpprobecombos.txt.gz";
////		String eqtlfile = "D:\\Work\\eqtlgeneremap\\test.txt";
//			String genomefasta = "D:\\Data\\HumanGenome\\human_g1k_v37.fasta.gz";
//			String ensemblannotation = "D:\\Data\\Ref\\Ensembl\\Homo_sapiens.GRCh37.71.gtf.gz";
//			String outdir = "D:\\Work\\eqtlgeneremap\\";
//			int windowsize = 2000000;
//			int readlen = 25;
//			int shift = 2;


        } else if (args[0].equals("align")) {

            if (args.length < 3) {
                System.out.println("Usage: align indir eqtlfile");
            } else {
                String indir = args[1];
                String eqtlfile = args[2];
                try {
                    mt.align(indir, eqtlfile);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        } else if (args[0].equals("quantifyIndividualEQTL")) {
            if (args.length < 5) {
                System.out.println("Usage quantifyIndividualEQTL eqtlfile indir gtfannotation out");
            } else {
                String eqtlfile = args[1];
                String indir = args[2];
                String ensemblannotation = args[3];
                String out = args[4];
                try {
                    mt.quantifyIndividualEQTL(eqtlfile, indir, ensemblannotation, out);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        } else if (args[0].equals("quantifysnp")) {
            if (args.length < 4) {
                System.out.println("Usage quantifysnp eqtlfile indir out");
            } else {
                String alignment = args[1];
                String genereadfile = args[2];
                String out = args[3];
                try {
                    mt.quantifySNP(alignment, genereadfile, out);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        } else if (args[0].equals("combine")) {

            if (args.length < 4) {
                System.out.println("Usage: combine alignout eqtl output");
            } else {

                String alignout = args[1];
                String eqtl = args[2];
                String output = args[3];
                try {
                    mt.combine(alignout, eqtl, output);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

        }



    }
}
