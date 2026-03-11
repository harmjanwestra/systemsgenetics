package transeqtlcrossmap;


import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;

import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.*;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.EQTL;
import umcg.genetica.io.trityper.QTLTextFile;
import umcg.genetica.text.Strings;

public class TransEQTLCrossMapBulk {

    // for future use
    public enum READCREATEMETHOD {
        FULLGENE,
        FIVEPRIME,
        THREEPRIME
    }

    public READCREATEMETHOD readcreatemethod;

    public void createSequences(String eqtlfile, String genomeFasta, String gtffile, String outputdir,
                                int ciswindowsize, int readlen, int shift, boolean exportindividualgenes,
                                boolean exportfullexonsequences, boolean overwriteexisitingsnps, String jobtemplate,
                                boolean useCisGeneAsCisWindow) throws IOException {

        System.out.println("--------------------------------");
        System.out.println("Cross mapping sequence generator");
        System.out.println("--------------------------------");
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        LocalDateTime now = LocalDateTime.now();
        System.out.println("Date: " + dtf.format(now));
        System.out.println("eQTL file: " + eqtlfile);
        System.out.println("Genome fasta: " + genomeFasta);
        System.out.println("GTF annotation: " + gtffile);
        System.out.println("Output dir: " + outputdir);
        System.out.println("Read length: " + readlen);
        System.out.println("Job template: " + jobtemplate);
        System.out.println("Use cis gene as reference window: " + useCisGeneAsCisWindow);
        System.out.println("(Additional) Cis window size: " + ciswindowsize);
        System.out.println();
        System.out.println();
        GTFAnnotation gtfAnnotation = new GTFAnnotation(gtffile);

        int batchsize = 20;
        Gpio.createDir(outputdir);
        TextFile logout = new TextFile(outputdir+"/makeseq.log",TextFile.W);
        String geneoutputdir = outputdir + "/testgenes/";
        String refwindowoutputdir = outputdir + "/refwindows/";
        Gpio.createDir(geneoutputdir);
        Gpio.createDir(refwindowoutputdir);

        System.out.println("Output of test genes: " + geneoutputdir);
        System.out.println("Output of reference windows: " + refwindowoutputdir);

        // create reference genome from snps
        System.out.println("Parsing QTL input: " + eqtlfile);
        TextFile eqtlfiletf = new TextFile(eqtlfile, TextFile.R);
        String[] header = eqtlfiletf.readLineElems(Strings.tab); // skip header
        int snpchrcol = -1;
        int snpposcol = -1;
        int snpnamecol = -1;
        int transgenenamecol = -1;
        int cisgenenamecol = -1;

        for (int i = 0; i < header.length; i++) {
            if (header[i].toLowerCase().equals("variantchr") || header[i].toLowerCase().equals("snpchr")) {
                snpchrcol = i;
            } else if (header[i].toLowerCase().equals("variantpos") || header[i].toLowerCase().equals("snppos")) {
                snpposcol = i;
            } else if (header[i].toLowerCase().equals("variantname") || header[i].toLowerCase().equals("snpname")) {
                snpnamecol = i;
            } else if (header[i].toLowerCase().equals("cisgene")) {
                cisgenenamecol = i;
            } else if (header[i].toLowerCase().equals("transgene") || header[i].toLowerCase().equals("cogene")) {
                transgenenamecol = i;
            }
        }

        System.out.println("Variant chromosome is on column: " + snpchrcol);
        System.out.println("Variant position is on column: " + snpposcol);
        System.out.println("Variant name is on column: " + snpnamecol);
        System.out.println("Cis gene ID is on column: " + cisgenenamecol);
        System.out.println("Trans gene ID is on column: " + transgenenamecol);

        if (snpchrcol == -1 || snpposcol == -1 || snpnamecol == -1 || transgenenamecol == -1) {
            System.out.println("Not all required columns are present in " + eqtlfile);
            System.out.println("Check above which ones are missing. Columns with the following names are required:");
            System.out.println("variantchr or snpchr");
            System.out.println("variantpos or snppos");
            System.out.println("variantname or snpname");
            System.out.println("cisgene");
            System.out.println("transgene or cogene");
            System.exit(-1);
        }

        ArrayList<Triple<String, Chromosome, Integer>> snps = new ArrayList<>();
        HashSet<String> snpids = new HashSet<>();
        HashSet<String> transgeneIds = new HashSet<>();
        HashSet<String> cisgeneIds = new HashSet<>();
        String[] elems = eqtlfiletf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            // gene_id variant chromosome      bp      eff_allele      non_eff_allele
            // ENSG00000 rs13    10      920        G       GAA
            // gene snpname  snpchr snppos
            Chromosome chr = Chromosome.parseChr(elems[snpchrcol]);
            Integer snppos = Integer.parseInt(elems[snpposcol]);
            String snpname = elems[snpnamecol].replaceAll(":", "_");
            String transgenename = elems[transgenenamecol];
            String cisgenename = elems[cisgenenamecol];
            if (!snpids.contains(snpname)) {
                Triple<String, Chromosome, Integer> snp = new Triple<>(snpname, chr, snppos);
                snps.add(snp);
                snpids.add(snpname);
            }

            transgeneIds.add(transgenename);
            cisgeneIds.add(cisgenename);
            elems = eqtlfiletf.readLineElems(TextFile.tab);
            ctr++;
//			if (ctr % 100 == 0) {
//				System.out.print("\r" + ctr + " lines parsed..");
//				break;
//			}
        }
        eqtlfiletf.close();

        System.out.println();
        System.out.println(snps.size() + " variants loaded");
        System.out.println(cisgeneIds.size() + " cis genes loaded");
        System.out.println(transgeneIds.size() + " trans/co-genes to export");

        logout.writeln("variants: "+snps.size()+", cis genes: "+cisgeneIds.size()+", trans genes: "+transgeneIds.size());

        exportTransGenes(genomeFasta, geneoutputdir, transgeneIds, gtfAnnotation, exportindividualgenes, exportfullexonsequences, readlen, shift, logout);
        exportRefSequences(snps, cisgeneIds, gtfAnnotation, genomeFasta, ciswindowsize, refwindowoutputdir, overwriteexisitingsnps, useCisGeneAsCisWindow, logout);
//
        // make index scripts
        TextFile tpl = new TextFile(jobtemplate, TextFile.R);
        ArrayList<String> templatelines = tpl.readAsArrayList();
        tpl.close();

        // create index slurm jobs
        String indexoutputdir = outputdir + "/1-indexjobs/";
        String indexlogoutputdir = outputdir + "/1-indexjobs/logs/";
        Gpio.createDir(indexoutputdir);
        Gpio.createDir(indexlogoutputdir);
        int bctr = 1;
        int cctr = 0;
        TextFile jobout = new TextFile(indexoutputdir + bctr + ".sh", TextFile.W);
        for (String line : templatelines) {
            if (line.contains("JOBNAME")) {
                jobout.writeln(line.replace("JOBNAME", "idx-" + bctr));
            } else if (line.contains("OUTLOGFILE")) {
                jobout.writeln(line.replace("OUTLOGFILE", indexlogoutputdir + "/idx-" + bctr + ".log"));
            } else if (line.contains("ERRLOGFILE")) {
                jobout.writeln(line.replace("ERRLOGFILE", indexlogoutputdir + "/idx-" + bctr + ".err"));
            } else if (!line.contains("CMD")) {
                jobout.writeln(line); // maximum JANK
            }
        }

        HashSet<String> idsToOutput = snpids;
        if (useCisGeneAsCisWindow) {
            idsToOutput = cisgeneIds;
        }
        for (String id : idsToOutput) {
            System.out.println(refwindowoutputdir + id + ".fa.gz");

            if (Gpio.exists(refwindowoutputdir + id + ".fa.gz")) {
                // "bwa index " + outputdir + "/snps/" + snp + ".fa.gz"
                for (String line : templatelines) {
                    if (line.contains("CMD")) {
                        String command = "bwa index " + refwindowoutputdir + id + ".fa.gz";
                        jobout.writeln(line.replaceAll("CMD", command));
                    }
                }
            }
            cctr++;
            if (cctr % batchsize == 0) {
                jobout.close();
                bctr += 1;
                jobout = new TextFile(indexoutputdir + bctr + ".sh", TextFile.W);
                for (String line : templatelines) {
                    if (line.contains("JOBNAME")) {
                        jobout.writeln(line.replace("JOBNAME", "idx-" + bctr));
                    } else if (line.contains("OUTLOGFILE")) {
                        jobout.writeln(line.replace("OUTLOGFILE", indexlogoutputdir + "/idx-" + bctr + ".log"));
                    } else if (line.contains("ERRLOGFILE")) {
                        jobout.writeln(line.replace("ERRLOGFILE", indexlogoutputdir + "/idx-" + bctr + ".err"));
                    } else if (!line.contains("CMD")) {
                        jobout.writeln(line);
                    }
                }
            }
        }

        jobout.close();

//
//        TextFile snplistout = new TextFile(refwindowoutputdir + "snps.txt", TextFile.W);
//        for (String snp : snpids) {
//            System.out.println(refwindowoutputdir + snp + ".fa.gz");
//            if (Gpio.exists(refwindowoutputdir + snp + ".fa.gz")) {
//                String ln = "bwa aln -t 16 -l 10 -0 ./snps/" + snp + ".fa.gz ./genes/allgenes.fa.gz > ./alignments/" + snp + "_allgenes.sai";
//                String ln2 = "bwa samse ./snps/" + snp + ".fa.gz ./alignments/" + snp + "_allgenes.sai ./genes/allgenes.fa.gz > ./alignments/" + snp + "_allgenes.sam";
//            }
//            snplistout.writeln(snp);
//        }
//        snplistout.close();


        // make alignment script
        eqtlfiletf.open();
        eqtlfiletf.readLine();
        elems = eqtlfiletf.readLineElems(TextFile.tab); // header

//        TextFile tfs2 = new TextFile(outputdir + "2-align.sh", TextFile.W);
//        TextFile tfs3 = new TextFile(outputdir + "3-samse.sh", TextFile.W);

        batchsize = 1000;
        String alignjobsoutputdir = outputdir + "/2-alignjobs/";
        String alignjobslogoutputdir = outputdir + "/2-alignjobs/logs/";
        Gpio.createDir(alignjobsoutputdir);
        Gpio.createDir(alignjobslogoutputdir);

        String samsejobsoutputdir = outputdir + "/3-samsejobs/";
        String samsejobslogoutputdir = outputdir + "/3-samsejobs/logs/";
        Gpio.createDir(samsejobsoutputdir);
        Gpio.createDir(samsejobslogoutputdir);

        String alignmentsdir = outputdir + "/alignments/";
        Gpio.createDir(alignmentsdir);


        bctr = 1;
        cctr = 0;
        TextFile alignjobout = new TextFile(alignjobsoutputdir + bctr + ".sh", TextFile.W);
        TextFile samsejobout = new TextFile(samsejobsoutputdir + bctr + ".sh", TextFile.W);
        for (String line : templatelines) {

            if (line.contains("JOBNAME")) {
                alignjobout.writeln(line.replace("JOBNAME", "aln-" + bctr));
                samsejobout.writeln(line.replace("JOBNAME", "sam-" + bctr));
            } else if (line.contains("OUTLOGFILE")) {
                alignjobout.writeln(line.replace("OUTLOGFILE", alignjobslogoutputdir + "/aln-" + bctr + ".log"));
                samsejobout.writeln(line.replace("OUTLOGFILE", samsejobslogoutputdir + "/sam-" + bctr + ".log"));
            } else if (line.contains("ERRLOGFILE")) {
                alignjobout.writeln(line.replace("ERRLOGFILE", alignjobslogoutputdir + "/aln-" + bctr + ".err"));
                samsejobout.writeln(line.replace("ERRLOGFILE", samsejobslogoutputdir + "/sam-" + bctr + ".err"));
            } else if (!line.contains("CMD")) {
                alignjobout.writeln(line);
                samsejobout.writeln(line);
            }
        }

        while (elems != null) {

            String snpname = elems[snpnamecol].replaceAll(":", "_");
            String cisgene = elems[cisgenenamecol];
            String transgene = elems[transgenenamecol];
            String refname = snpname;
            if (useCisGeneAsCisWindow) {
                refname = cisgene;
            }
            if (Gpio.exists(refwindowoutputdir + refname + ".fa.gz") && Gpio.exists(geneoutputdir + transgene + ".fa.gz")) {

//                TextFile alignjobout = new TextFile(alignjobsoutputdir + snp + "_" + gene + ".sh", TextFile.W);
                for (String line : templatelines) {
                    if (line.contains("CMD")) {
                        String command = "bwa aln -t 1 -l 10 -0 " + refwindowoutputdir + refname + ".fa.gz " + geneoutputdir + transgene + ".fa.gz > " + alignmentsdir + refname + "_" + transgene + ".sai";
                        String outln = line.replaceAll("CMD", command);
                        alignjobout.writeln(outln);
                        command = "bwa samse " + refwindowoutputdir + refname + ".fa.gz " + alignmentsdir + refname + "_" + transgene + ".sai " + geneoutputdir + transgene + ".fa.gz > " + alignmentsdir + refname + "_" + transgene + ".sam";
                        outln = line.replaceAll("CMD", command);
                        samsejobout.writeln(outln);
                    }
                }

                cctr++;
                if (cctr % batchsize == 0) {
                    bctr++;
                    alignjobout.close();
                    samsejobout.close();

                    alignjobout = new TextFile(alignjobsoutputdir + bctr + ".sh", TextFile.W);
                    samsejobout = new TextFile(samsejobsoutputdir + bctr + ".sh", TextFile.W);


                    for (String line : templatelines) {
                        if (line.contains("JOBNAME")) {
                            alignjobout.writeln(line.replace("JOBNAME", "aln-" + bctr));
                            samsejobout.writeln(line.replace("JOBNAME", "sam-" + bctr));
                        } else if (line.contains("OUTLOGFILE")) {
                            alignjobout.writeln(line.replace("OUTLOGFILE", alignjobslogoutputdir + "/aln-" + bctr + ".log"));
                            samsejobout.writeln(line.replace("OUTLOGFILE", samsejobslogoutputdir + "/sam-" + bctr + ".log"));
                        } else if (line.contains("ERRLOGFILE")) {
                            alignjobout.writeln(line.replace("ERRLOGFILE", alignjobslogoutputdir + "/aln-" + bctr + ".err"));
                            samsejobout.writeln(line.replace("ERRLOGFILE", samsejobslogoutputdir + "/sam-" + bctr + ".err"));
                        } else if (!line.contains("CMD")) {
                            alignjobout.writeln(line);
                            samsejobout.writeln(line);
                        }
                    }

                }

            }
            elems = eqtlfiletf.readLineElems(TextFile.tab);
        }
        eqtlfiletf.close();
        alignjobout.close();
        samsejobout.close();
//        tfs2.close();
//        tfs3.close();
        logout.close();
    }

    private void exportRefSequences(ArrayList<Triple<String, Chromosome, Integer>> snps,
                                    HashSet<String> cisgeneIds, GTFAnnotation gtfAnnotation, String genomeFasta,
                                    int ciswindowsize, String refwindowoutputdir, boolean overwriteexisitingsnps,
                                    boolean useCisGeneAsCisWindow, TextFile logout) {
        System.out.println("Loading genome from: " + genomeFasta);
        FastaSequenceFile fastaFile = new FastaSequenceFile(new File(genomeFasta), false);

        ReferenceSequence seq = fastaFile.nextSequence();


        ExecutorService ex = Executors.newWorkStealingPool(16);
        AtomicInteger nrsnpsexported = new AtomicInteger();

        ArrayList<Future> taskoutput = new ArrayList<>();
        while (seq != null) {
            String seqname = seq.getName();
            while (seqname.contains("  ")) {
                seqname = seqname.replaceAll("  ", " ");
            }
            String[] seqnamelems = seqname.split(" ");
            seqname = seqnamelems[0];


            Chromosome seqchr = Chromosome.parseChr(seqname);

            if (seqchr.isAutosome()) {
                System.out.println("Getting sequence: " + seq.getName() + "\t" + seq.length() + "\tParsed seq: " + seqname);
                Runnable t = new ReferenceExportTask(seq.getBases(), Chromosome.parseChr(seqname), snps, cisgeneIds, gtfAnnotation, ciswindowsize, nrsnpsexported, refwindowoutputdir, overwriteexisitingsnps, useCisGeneAsCisWindow);

                Future<?> output = ex.submit(t);
                taskoutput.add(output);
            }
            seq = fastaFile.nextSequence();
        }


        boolean alldone = false;
        while (!alldone) {
            alldone = true;
            for (Future f : taskoutput) {
                if (!f.isDone()) {
                    alldone = false;
                }
            }
        }

        fastaFile.close();
        ex.shutdown();
    }

    private String byteToStr(byte[] bytes) throws UnsupportedEncodingException {
        return new String(bytes, StandardCharsets.UTF_8);
    }

    class ReferenceExportTask implements Runnable {
        private final boolean overwriteexisitingsnps;
        private final byte[] bases;
        private final Chromosome seqchr;
        private final ArrayList<Triple<String, Chromosome, Integer>> snps;
        private final int windowsize;
        private final AtomicInteger nrsnpsexported;
        private final String snpoutputdir;
        private final GTFAnnotation gtf;
        private final HashSet<String> geneIDs;
        private final boolean useCisGeneAsCisWindow;

        public ReferenceExportTask(byte[] bases,
                                   Chromosome seqchr,
                                   ArrayList<Triple<String, Chromosome, Integer>> snps,
                                   HashSet<String> geneIDs,
                                   GTFAnnotation gtf,
                                   int windowsize,
                                   AtomicInteger nrsnpsexported,
                                   String snpoutputdir,
                                   boolean overwriteexisitingsnps,
                                   boolean useCisGeneAsCisWindow) {
            this.bases = bases;
            this.gtf = gtf;
            this.geneIDs = geneIDs;
            this.seqchr = seqchr;
            this.snps = snps;
            this.windowsize = windowsize;
            this.nrsnpsexported = nrsnpsexported;
            this.snpoutputdir = snpoutputdir;
            this.overwriteexisitingsnps = overwriteexisitingsnps;
            this.useCisGeneAsCisWindow = useCisGeneAsCisWindow;
        }

        @Override
        public void run() {
            try {

                // determine what to export
                ArrayList<Feature> features = new ArrayList<>();
                int nrToExport = snps.size();
                if (useCisGeneAsCisWindow) {
                    for (String geneStr : geneIDs) {
                        Gene geneObj = gtf.getStrToGene().get(geneStr);
                        if (geneObj.getChromosome().equals(seqchr)) {
                            features.add(geneObj);
                        }
                    }
                    nrToExport = geneIDs.size();
                } else {
                    for (Triple<String, Chromosome, Integer> snp : snps) {
                        if (snp.getMiddle().equals(seqchr)) {
                            Feature f = new SNPFeature(seqchr, snp.getRight(), snp.getRight());
                            f.setName(snp.getLeft());
                            features.add(f);
                        }
                    }
                }

                for (Feature f : features) {

                    String filename = snpoutputdir + f.getName() + ".fa.gz";
                    if (Gpio.exists(filename) && !overwriteexisitingsnps && Gpio.getFileSize(filename) > 0) {
                        // do not overwrite
                        nrsnpsexported.getAndIncrement();
                    } else {
                        int posLeft = f.getStart();
                        int posRight = f.getStop();

                        int windowleft = posLeft - windowsize;
                        int windowright = posRight + windowsize;

                        if (windowleft < 0) {
                            windowleft = 0;
                        }
                        if (windowright > bases.length - 1) {
                            windowright = bases.length - 1;
                        }

                        int nrbases = windowright - windowleft;
                        if (nrbases <= 0) {
                            System.err.println("Error: <=0 bases..");
                            System.err.println(nrsnpsexported.get() + "/" + snps.size() + " | Ref: " + f.getName() + " " + f.toString() + " left: " + windowleft + " right: " + windowright + " size: " + nrbases);
                            System.exit(-1);
                        }

                        byte[] basewindow = new byte[nrbases];
                        System.arraycopy(bases, windowleft, basewindow, 0, nrbases);

                        // make into fasta
                        TextFile snpfastaout = new TextFile(filename, TextFile.W);
                        String fastaStr = null;
                        try {
                            int featlen = f.getSize();
                            fastaStr = ">" + f.getName() + "_" + f.toString() + "|featureLen:" + featlen + "|nrBases:" + nrbases
                                    + "\n" + byteToStr(basewindow);
                        } catch (UnsupportedEncodingException e) {
                            e.printStackTrace();
                        }
                        snpfastaout.writeln(fastaStr);
                        snpfastaout.close();
                        System.out.println(nrsnpsexported + "/" + nrToExport + " | Ref: " + f.getName() + " " + f.toString() + " left: " + windowleft + " right: " + windowright + " size: " + nrbases);
                        nrsnpsexported.getAndIncrement();
                    }
                }
            } catch (
                    Exception e) {
                e.printStackTrace();
            }
        }
    }

    protected void exportTransGenes(String genomeFasta, String geneoutputdir,
                                    HashSet<String> geneIds, GTFAnnotation gtfAnnotation,
                                    boolean exportindividualgenes, boolean exportfullexonsequences,
                                    int readlen, int shift, TextFile logout) throws IOException {


        // convert gene strings to ensembl gene ids
        HashMap<String, Gene> geneHash = gtfAnnotation.getStrToGene();
        ArrayList<Gene> genesToExport = new ArrayList<>();
        for (String g : geneIds) {
            Gene gobj = geneHash.get(g);
            if (gobj != null) {
                genesToExport.add(gobj);
            } else {
                logout.writeln(g+" not found in GTF annotation.");
            }
        }

        System.out.println(genesToExport.size() + " genes found in annotation.");
        logout.writeln(genesToExport.size() + " genes found in annotation.");

        if (genesToExport.isEmpty()) {
            System.err.println("No genes to export..?");
            System.exit(-1);
        }

        // create 1mb regions around snps
        System.out.println("Loading genome from: " + genomeFasta);
        FastaSequenceFile fastaFile = new FastaSequenceFile(new File(genomeFasta), false);
        ReferenceSequence seq = fastaFile.nextSequence();
        int nrsnpsexported = 0;
        int nrgenesexported = 0;
        TextFile allgenefa = new TextFile(geneoutputdir + "allgenes.fa.gz", TextFile.W);
        TextFile allgenefaignoringexons = new TextFile(geneoutputdir + "allgenes_ignoringexons.fa.gz", TextFile.W);
        TextFile geneStats = new TextFile(geneoutputdir + "allgenes.stats.txt", TextFile.W);
        geneStats.writeln("Gene\tnrOfReads\tnrOfBases\tnrTranscripts\tnrExons\tnrOfReadsIgnoringExons\tnrOfBasesIgnoringExons");
        while (seq != null) {

            byte[] bases = seq.getBases();

            String seqname = seq.getName();
            while (seqname.contains("  ")) {
                seqname = seqname.replaceAll("  ", " ");
            }
            String[] seqnamelems = seqname.split(" ");
            seqname = seqnamelems[0];
            System.out.println("Getting sequence: " + seq.getName() + "\t" + seq.length() + "\tParsed seq: " + seqname);

            Chromosome seqchr = Chromosome.parseChr(seqname);


            System.out.println("Now exporting genes");
            // now parse the ensemblgenes on this chr
            for (Gene g : genesToExport) {
                if (g.getChromosome().equals(seqchr)) {
                    TextFile genesfastaout = null;
                    if (exportindividualgenes) {
                        genesfastaout = new TextFile(geneoutputdir + g.getName() + ".fa.gz", TextFile.W);
                    }
                    ArrayList<Transcript> transcripts = g.getTranscripts();
                    HashSet<Exon> allexons = new HashSet<Exon>();
                    for (Transcript t : transcripts) {
                        ArrayList<Exon> exons = t.getExons();
                        for (Exon e : exons) {
                            if (e.getType().equals(FeatureType.EXON)) {
                                allexons.add(e);
                            }
                        }
                    }
//					System.out.println("Gene: " + g.getName() + "\tTotal exons: " + allexons.size());


                    if (exportfullexonsequences) {
                        TextFile genesfastaoutFull = new TextFile(geneoutputdir + g.getName() + "_exons.fa.gz", TextFile.W);
                        for (Exon currentexon : allexons) {
                            if (currentexon.getType().equals(FeatureType.EXON)) {
                                int len = currentexon.getStop() - currentexon.getStart();
                                byte[] b = new byte[len];
                                System.arraycopy(bases, currentexon.getStart(), b, 0, len);
                                String fastaStr = ">" + g.getName() + "_" + currentexon.getName() + "_" + currentexon.getChromosome().toString() + "_" + currentexon.getStart() + "_" + currentexon.getStop()
                                        + "\n" + byteToStr(b);
                                genesfastaoutFull.writeln(fastaStr);
                            }
                        }
                        genesfastaoutFull.close();
                    }

                    // split up the transcript in readlen reads
                    int nrseqs = 0;
                    for (Transcript t : transcripts) {
                        ArrayList<Exon> exons = t.getExons();

                        // get only exons
                        {
                            ArrayList<Exon> tmpexons = new ArrayList<>();
                            for (Exon e : exons) {
                                if (e.getType().equals(FeatureType.EXON)) {
                                    tmpexons.add(e);
                                }
                            }
                            exons = tmpexons;
                        }

                        Collections.sort(exons, new FeatureComparator(true)); // sort the exons

//						if (t.getStrand().equals(Strand.NEG)) {
//							System.out.println("Got one...");
//							for (int i = 0; i < exons.size(); i++) {
//								System.out.println(i + "\t" + exons.get(i).toString());
//							}
//							System.exit(-1);
//						}

                        // output exon sequences
                        int exonctr = 0;
                        int tstart = exons.get(0).getStart();

//						if (exons.size() > 1) {
//							System.out.println("Bench this");
//						}

                        int tstop = exons.get(exons.size() - 1).getStop();
                        int currentPos = tstart;
                        Exon currentexon = exons.get(0);
                        while (currentPos < tstop) {

                            // iterate current exon
                            while (currentPos + readlen <= currentexon.getStop()) {
                                byte[] b = new byte[readlen];
                                System.arraycopy(bases, currentPos, b, 0, readlen);
                                String fastaStr = ">" + g.getName() + "_" + t.getName() + "_" + currentexon.getName() + "_" + currentexon.getChromosome().toString() + "_" + currentPos + "_" + (currentPos + readlen)
                                        + "\n" + byteToStr(b);
                                if (exportindividualgenes) {
                                    genesfastaout.writeln(fastaStr);
                                }
                                allgenefa.writeln(fastaStr);
// System.out.println(nrseqs + "\texon: " + exonctr + "\tsta: " + currentPos + "\tsto: " + (currentPos + readlen) + "\texonstop: " + currentexon.getStop());
                                nrseqs++;
                                currentPos += shift;

                            }

                            // process exon exon boundary

                            // if there are exons remaining
                            if (exonctr + 1 < exons.size()) {
                                while (currentPos < currentexon.getStop()) { // iterate until we combineEPRSWithCisAndTransFX out of currentexon
                                    int basesFromCurrentExon = currentexon.getStop() - currentPos;
                                    byte[] b = new byte[readlen];
                                    System.arraycopy(bases, currentPos, b, 0, basesFromCurrentExon);

                                    int basesRemaining = readlen - basesFromCurrentExon;
                                    // now we should check whether the next exon actually has enough bases
                                    int tmpexonctr = exonctr + 1;
                                    while (basesRemaining > 0 && tmpexonctr < exons.size()) {
                                        // get remaining bases from next exons
                                        Exon tmpexon = exons.get(tmpexonctr);
                                        int tmpexonstart = tmpexon.getStart();
                                        int tmpexonstop = tmpexon.getStop();

                                        int basesToTakeFromExon = basesRemaining;
                                        if (tmpexonstart + basesRemaining > tmpexonstop) {
                                            basesToTakeFromExon = tmpexonstop - tmpexonstart;
                                        }


                                        System.arraycopy(bases, tmpexonstart, b, basesFromCurrentExon, basesToTakeFromExon);
                                        basesRemaining -= basesToTakeFromExon;
//										System.out.println("Boundary: " + nrseqs + "\texon: " + exonctr + "\ttmpexon: " + tmpexonctr + "\tsta: " + currentPos + "\texonstop: " + currentexon.getStop() + "\t" + basesRemaining);
                                        tmpexonctr++;
                                    }
                                    if (basesRemaining == 0) {
                                        // write the resulting freakshow of a sequence
                                        String fastaStr = ">" + g.getName() + "_" + t.getName() + "_" + currentexon.getName() + "_" + currentPos + "_" + (currentPos + readlen)
                                                + "\n" + byteToStr(b);
                                        if (exportindividualgenes) {
                                            genesfastaout.writeln(fastaStr);
                                        }
                                        allgenefa.writeln(fastaStr);
                                        nrseqs++;
                                    }

                                    currentPos += shift;
                                }
                            }

                            // start at next exon
                            if (exonctr < exons.size()) {
//								System.out.println("Investigating exon: " + exonctr);
                                currentexon = exons.get(exonctr);
                                currentPos = currentexon.getStart();
                            } else {
//								System.out.println("Reached end of transcript");
                                currentPos = currentexon.getStop();
                            }
                            exonctr++;
                        }

                    }

                    if (exportindividualgenes) {
                        genesfastaout.close();
                    }

                    // write the sequence, while ignoring exons
                    int sta = g.getStart();
                    int sto = g.getStop();

                    int currentPos = sta;
                    while (currentPos < sto) {
                        byte[] b = new byte[readlen];
                        System.arraycopy(bases, currentPos, b, 0, readlen);
                        String fastaStr = ">" + g.getName() + "\t" + g.getChromosome().toString() + "_" + currentPos + "_" + (currentPos + readlen)
                                + "\n" + byteToStr(b);
                        allgenefaignoringexons.writeln(fastaStr);
                        currentPos += shift;
                    }


                    String genestatout = g.getName() + "\t" + nrseqs + "\t" + (nrseqs * readlen) + "\t" + transcripts.size() + "\t" + allexons.size();
                    geneStats.writeln(genestatout);

                    System.out.println(nrgenesexported + "/" + genesToExport.size() + " | Gene: " + g.getName() + " Chr: " + g.getChromosome().toString() + " start: " + g.getStart() + " stop: " + g.getStop() + " Strand: " + g.getStrand() + " nr sequences: " + nrseqs);
                    nrgenesexported++;
                }
            }
            System.out.println("Done exporting genes");
            seq = fastaFile.nextSequence();
        }
        geneStats.close();
        allgenefa.close();
    }

    public void combine(String alignmentoutput, String eqtlfile, String output) throws IOException {

        HashMap<String, String> snpprobeToOutput = new HashMap<String, String>();
        TextFile tf = new TextFile(alignmentoutput, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        while (elems != null) {
            snpprobeToOutput.put(elems[0] + "_" + elems[1], Strings.concat(elems, Strings.tab, 2, elems.length));
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile out = new TextFile(output, TextFile.W);
        QTLTextFile qtf = new QTLTextFile(eqtlfile, QTLTextFile.R);
        out.writeln("SNP\tSNPChr\tSNPChrPos\tGene\tGeneChr\tGeneChrPos\tZ\tAbsZ\tNrReads\t%BasesMapped\t%ReadsMapped");

        Iterator<EQTL> it = qtf.getEQtlIterator();
        while (it.hasNext()) {
            EQTL e = it.next();
            String query = e.getRsName() + "_" + e.getProbe();
            String append = snpprobeToOutput.get(query);

            String lnout = e.getRsName()
                    + "\t" + e.getRsChr()
                    + "\t" + e.getRsChrPos()
                    + "\t" + e.getProbe()
                    + "\t" + e.getProbeChr()
                    + "\t" + e.getProbeChrPos()
                    + "\t" + e.getZscore()
                    + "\t" + e.getZscoreAbs()
                    + "\t" + append;
            out.writeln(lnout);
        }
        tf.close();
        out.close();

    }

    public void quantifyIndividualEQTL(String eqtlfile, String indir, String ensemblannotation,
                                       String out, boolean useCisGeneAsCisWindow, int ciswindowsize) throws IOException {

        System.out.println("----------------------------");
        System.out.println("Cross mapping quantification");
        System.out.println("----------------------------");
        DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        LocalDateTime now = LocalDateTime.now();

        System.out.println("QTL file: " + eqtlfile);
        System.out.println("Indir: " + indir);
        System.out.println("GTF annotation: " + ensemblannotation);
        System.out.println("Out: " + out);
        System.out.println("Use cis gene as reference window: " + useCisGeneAsCisWindow);
        System.out.println("(Additional) Cis window size: " + ciswindowsize);
        System.out.println();
        System.out.println();

        GTFAnnotation s = new GTFAnnotation(ensemblannotation);

        // create reference genome from snps
        System.out.println("Parsing QTL input: " + eqtlfile);
        TextFile eqtlfiletf = new TextFile(eqtlfile, TextFile.R);
        String[] header = eqtlfiletf.readLineElems(Strings.tab); // skip header
        int snpchrcol = -1;
        int snpposcol = -1;
        int snpnamecol = -1;
        int transgenenamecol = -1;
        int cisgenenamecol = -1;

        for (int i = 0; i < header.length; i++) {
            if (header[i].toLowerCase().equals("variantchr") || header[i].toLowerCase().equals("snpchr")) {
                snpchrcol = i;
            } else if (header[i].toLowerCase().equals("variantpos") || header[i].toLowerCase().equals("snppos")) {
                snpposcol = i;
            } else if (header[i].toLowerCase().equals("variantname") || header[i].toLowerCase().equals("snpname")) {
                snpnamecol = i;
            } else if (header[i].toLowerCase().equals("cisgene")) {
                cisgenenamecol = i;
            } else if (header[i].toLowerCase().equals("transgene") || header[i].toLowerCase().equals("cogene")) {
                transgenenamecol = i;
            }
        }

        System.out.println("Variant chromosome is on column: " + snpchrcol);
        System.out.println("Variant position is on column: " + snpposcol);
        System.out.println("Variant name is on column: " + snpnamecol);
        System.out.println("Cis gene ID is on column: " + cisgenenamecol);
        System.out.println("Trans gene ID is on column: " + transgenenamecol);

        if (snpchrcol == -1 || snpposcol == -1 || snpnamecol == -1 || transgenenamecol == -1) {
            System.out.println("Not all required columns are present in " + eqtlfile);
            System.out.println("Check above which ones are missing. Columns with the following names are required:");
            System.out.println("variantchr or snpchr");
            System.out.println("variantpos or snppos");
            System.out.println("variantname or snpname");
            System.out.println("cisgene");
            System.out.println("transgene or cogene");
            System.exit(-1);
        }

        ArrayList<Triple<String, Chromosome, Integer>> snps = new ArrayList<>();
        HashSet<String> snpids = new HashSet<>();
        HashSet<String> transgeneIds = new HashSet<>();
        HashSet<String> cisgeneIds = new HashSet<>();
        String[] elems = eqtlfiletf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            // gene_id variant chromosome      bp      eff_allele      non_eff_allele
            // ENSG00000 rs13    10      920        G       GAA
            // gene snpname  snpchr snppos
            Chromosome chr = Chromosome.parseChr(elems[snpchrcol]);
            Integer snppos = Integer.parseInt(elems[snpposcol]);
            String snpname = elems[snpnamecol].replaceAll(":", "_");
            String transgenename = elems[transgenenamecol];
            String cisgenename = elems[cisgenenamecol];
            if (!snpids.contains(snpname)) {
                Triple<String, Chromosome, Integer> snp = new Triple<>(snpname, chr, snppos);
                snps.add(snp);
                snpids.add(snpname);
            }

            transgeneIds.add(transgenename);
            cisgeneIds.add(cisgenename);
            elems = eqtlfiletf.readLineElems(TextFile.tab);
            ctr++;
//			if (ctr % 100 == 0) {
//				System.out.print("\r" + ctr + " lines parsed..");
//				break;
//			}
        }
        eqtlfiletf.close();

        System.out.println("Parsing gene fasta files:");
        HashMap<String, Integer> geneToNrReads = new HashMap<String, Integer>();
        int gctr = 0;
        int missing = 0;
        for (String genename : transgeneIds){
            String fasta = indir + "/testgenes/" + genename + ".fa.gz";
            if(Gpio.exists(fasta)) {
                TextFile faIN = new TextFile(fasta, TextFile.R);
                int nrlines = faIN.countLines();
                nrlines /= 2;
                geneToNrReads.put(genename, nrlines);
                gctr++;
                if (gctr % 1000 == 0) {
                    System.out.print(gctr + " genes parsed\r");
                }
            } else {
                System.out.println("Warning: cannot find "+fasta);
                missing++;
            }
        }
        System.out.println(gctr + " genes parsed");
        System.out.println(missing +" genes missing.");


        System.out.println("Parsing alignments...");
        // iterate the QTL file once more
        eqtlfiletf.open();
        eqtlfiletf.readLine(); // skip header
        elems = eqtlfiletf.readLineElems(TextFile.tab);
        ctr = 0;

        TextFile outtf = new TextFile(out, TextFile.W);

        String closeStr = "DistanceBetweenSNPAndClosestAlignment";
        if (useCisGeneAsCisWindow) {
            closeStr = "LeftMostAlignmentRelativeToWindow";
        }
        outtf.writeln("CisChr\t" +
                "BaseWindowStart\t" +
                "BaseWindowEnd\t" +
                "ReferenceSize\t" +
                "CisGene\t" +
                "CisGeneSymbol\t" +
                "CisGeneStart\t" +
                "CisGeneEnd\t" +
                "CisGeneStrand\t" +
                "Variant\t" +
                "VariantPos\t" +
                "TransGene\t" +
                "TransGeneSymbol\t" +
                "TransGeneChr\t" +
                "TransGeneStart\t" +
                "TransGeneEnd\t" +
                "TransGeneStrand\t" +
                "TransGeneSameChrAsSNPorCisGene\t" +
                "NrReads\t" +
                "NumberOfReadBasesMapped\t" +
                "ProportionReadsMapped\t"
                + closeStr);

        int missingalignments = 0;
        int parsedalignments = 0;
        while (elems != null) {
            Chromosome snpchr = Chromosome.parseChr(elems[snpchrcol]);
            Integer snppos = Integer.parseInt(elems[snpposcol]);
            String snpname = elems[snpnamecol].replaceAll(":", "_");
            String transgenename = elems[transgenenamecol];
            String cisgenename = elems[cisgenenamecol];

            Gene transGeneObj = s.getStrToGene().get(transgenename);
            Gene cisGeneObj = s.getStrToGene().get(cisgenename);

            // get the alignment
            String alignment = indir + "/alignments/" + snpname + "_" + transgenename + ".sam";
            if (useCisGeneAsCisWindow) {
                alignment = indir + "/alignments/" + cisgenename + "_" + transgenename + ".sam";
            }

            if (Gpio.exists(alignment)) {
                parsedalignments++;
//				System.out.println(alignment);
                SAMFileReader reader = new SAMFileReader(new File(alignment));

                SAMSequenceDictionary seqDict = reader.getFileHeader().getSequenceDictionary();
                long seqlen = seqDict.getReferenceLength();
                int relativesnppos = (int) seqlen / 2;
                int closest = Integer.MAX_VALUE;
                SAMRecordIterator it = reader.iterator();
                double totalperc = 0;


                // iterate all alignments, take the ones that are properly mapping
                while (it.hasNext()) {
                    net.sf.samtools.SAMRecord record = it.next();
                    List<AlignmentBlock> blocks = record.getAlignmentBlocks();
                    int readlen = record.getReadLength();
                    int mappedlen = 0;

                    int start = record.getAlignmentStart();
                    int relpos = Math.abs(relativesnppos - start);
                    if (relpos < closest) {
                        closest = relpos;
                    }


                    for (AlignmentBlock b : blocks) {
                        mappedlen += b.getLength();

                    }

                    double perc = (double) mappedlen / readlen;
                    totalperc += perc;
                }

                int nrReadsTotal = geneToNrReads.get(transgenename);
                // System.out.println(alignment + "\t" + relativesnppos + "\t" + closest);
                double propreads = (totalperc / nrReadsTotal);
                if (Double.isNaN(propreads)) {
                    propreads = -1;
                }


                int windowstart = snppos - ciswindowsize;
                int windowstop = snppos + ciswindowsize;
                if (useCisGeneAsCisWindow) {
                    windowstart = cisGeneObj.getStart() - ciswindowsize;
                    windowstop = cisGeneObj.getStop() + ciswindowsize;
                }

                if(windowstart < 0){
                    windowstart = 0;
                }


                outtf.writeln(
                        snpchr
                                + "\t" + windowstart
                                + "\t" + windowstop
                                + "\t" + seqlen
                                + "\t" + cisgenename
                                + "\t" + cisGeneObj.getGeneSymbol()
                                + "\t" + cisGeneObj.getStart()
                                + "\t" + cisGeneObj.getStop()
                                + "\t" + cisGeneObj.getStrand()
                                + "\t" + snpname
                                + "\t" + snppos
                                + "\t" + transgenename
                                + "\t" + transGeneObj.getGeneSymbol()
                                + "\t" + transGeneObj.getChromosome()
                                + "\t" + transGeneObj.getStart()
                                + "\t" + transGeneObj.getStop()
                                + "\t" + transGeneObj.getStrand()
                                + "\t" + (transGeneObj.getChromosome().equals(snpchr))
                                + "\t" + nrReadsTotal
                                + "\t" + totalperc
                                + "\t" + propreads
                                + "\t" + closest);

                reader.close();
            } else {
                System.out.println("Could not find alignment: " + alignment);
                missingalignments++;
            }


            elems = eqtlfiletf.readLineElems(TextFile.tab);
            ctr++;
            if (ctr % 100 == 0) {
                System.out.print("\r" + ctr + " lines parsed..");
            }
        }
        eqtlfiletf.close();
        outtf.close();
        System.out.println(missingalignments + " missing alignments");
        System.out.println(parsedalignments +" parsed alignments");

    }

    public void quantifySNP(String alignmentfile, String geneReadFile, String out) throws IOException {


        // map gene name to nr of reads
        HashMap<String, Double> geneToReads = new HashMap<String, Double>();
        HashMap<String, Double> geneToBases = new HashMap<String, Double>();
        HashMap<String, Integer> geneIndex = new HashMap<String, Integer>();
        ArrayList<String> genenames = new ArrayList<>();
        TextFile tf = new TextFile(geneReadFile, TextFile.R);
        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;

        while (elems != null) {

            String gene = elems[0];
            Double nrOfReads = Double.parseDouble(elems[1]);
            Double nrOfBases = Double.parseDouble(elems[2]);

            if (!geneIndex.containsKey(gene)) {
                geneToReads.put(gene, nrOfReads);
                geneToBases.put(gene, nrOfBases);
                geneIndex.put(gene, ctr);
                genenames.add(gene);
                ctr++;
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        TextFile outtf = new TextFile(out, TextFile.W);
        String header = "Gene\tNrOfReadsTotal\tNrOfReadsMapped\tNrOfBasesTotal\tNrOfBasesMapped";
        outtf.writeln(header);

        // iterate the list of snps

        if (!Gpio.exists(alignmentfile)) {
            System.out.println("Expected file: " + alignmentfile + " but wasn't present on disk.");
        } else {


            double[] nrReadsMappedPerGene = new double[geneIndex.size()];
            double[] nrBasesMappedPerGene = new double[geneIndex.size()];

            SAMFileReader reader = new SAMFileReader(new File(alignmentfile));
            SAMRecordIterator it = reader.iterator();

            while (it.hasNext()) {
                net.sf.samtools.SAMRecord record = it.next();
                String recordname = record.getReadName();
                String gene = recordname.split("_")[0];
                Integer index = geneIndex.get(gene);
                if (index != null) {
                    List<AlignmentBlock> blocks = record.getAlignmentBlocks();
                    if (!blocks.isEmpty()) {
                        nrReadsMappedPerGene[index]++;
                        int mappedlen = 0;
                        for (AlignmentBlock b : blocks) {
                            mappedlen += b.getLength();
                        }
                        nrBasesMappedPerGene[index] += mappedlen;
                    }
                }
            }

            for (int g = 0; g < nrBasesMappedPerGene.length; g++) {
                String gene = genenames.get(g);
                double nrOfReads = geneToReads.get(gene);
                double nrOfBases = geneToBases.get(gene);
                double percentageOfReadsMapped = nrOfReads / nrReadsMappedPerGene[g];
                double percentageOfBasesMapped = nrOfBases / nrBasesMappedPerGene[g];
//				String outln = genenames.get(g) + "\t" + nrOfReads + "\t" + nrReadsMappedPerGene[g] + "\t" + percentageOfReadsMapped + "\t" + nrOfBasesTotal + "\t" + nrBasesMappedPerGene[g] + "\t" + percentageOfBasesMapped;
            }

            System.exit(-1);
//			outtf.writeln(outln);

            it.close();
            reader.close();

        }


        outtf.close();
    }


    public void align(String indir, String eqtlfile) throws IOException {

        ExecutorService ex = Executors.newFixedThreadPool(16);

        TextFile tf = new TextFile(eqtlfile, TextFile.R);

        tf.readLine();
        int nrAlignments = tf.countLines();
        tf.close();
        tf.open();

        tf.readLine();
        String[] elems = tf.readLineElems(TextFile.tab);
        int nrsubmitted = 0;
        AtomicInteger nrdone = new AtomicInteger();
        ProgressBar pb = new ProgressBar(nrAlignments);
        while (elems != null) {

            String snpname = elems[0];
            Gpio.createDir(indir + "/alignments/" + snpname + "/");
            String genename = elems[3];

            AlignmentTask t = new AlignmentTask();
            t.gene = genename;
            t.snp = snpname;
            t.dir = indir;
            t.nrdone = nrdone;

            ex.submit(t);
            nrsubmitted++;

            // don't let the buffer combineEPRSWithCisAndTransFX too full
//			if (nrdone.get() - nrsubmitted > 128) {
//				try {
//					Thread.sleep(1000);
//				} catch (InterruptedException e) {
//					e.printStackTrace();+
//				}
//			}

//			if (nrsubmitted == 10000) {
//				while (nrdone.get() < nrsubmitted) {
//					try {
//
//						pb.set(nrdone.get());
//						Thread.sleep(1000);
//					} catch (InterruptedException e) {
//						e.printStackTrace();
//					}
//				}
//				System.exit(-1);
//			}

            if (nrsubmitted % 1000 == 0) {
                pb.set(nrdone.get());
            }

            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

        while (nrdone.get() < nrsubmitted) {
            try {
                pb.set(nrdone.get());
//				System.out.print("\r" + nrdone.get() + "/" + nrsubmitted);
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        ex.shutdown();
        pb.close();

    }


    public class AlignmentTask implements Runnable {

        String snp;
        String gene;
        String dir;
        public AtomicInteger nrdone;

        @Override
        public void run() {


//			String ln = "bwa aln -l 10 -0 " + dir + "/snps/" + snp + ".fa.gz " + dir + "/genes/" + gene + ".fa.gz > " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai";
            String[] cmd1 = new String[]{
                    "bwa",
                    "aln",
                    "-l 10",
                    "-0",
                    dir + "/snps/" + snp + ".fa.gz",
                    dir + "/genes/" + gene + ".fa.gz"

            };


//			String ln2 = "bwa samse " + dir + "/snps/" + snp + ".fa.gz " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai " + dir + "/genes/" + gene + ".fa.gz > " + dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sam";
            String[] cmd2 = new String[]{
                    "bwa",
                    "samse",
                    dir + "/snps/" + snp + ".fa.gz",
                    dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai",
                    dir + "/genes/" + gene + ".fa.gz"

            };
            try {
                ProcessBuilder pb = new ProcessBuilder(cmd1);
                pb.redirectOutput(new File(dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sai"));
                Process process = pb.start();

                int errCode = process.waitFor();
//				System.out.println("Echo command executed, any errors? " + (errCode == 0 ? "No" : "Yes"));
                ProcessBuilder pb2 = new ProcessBuilder(cmd2);
                pb2.redirectOutput(new File(dir + "/alignments/" + snp + "/" + snp + "_" + gene + ".sam"));
                Process process2 = pb2.start();
                int errCode2 = process2.waitFor();
//				System.out.println("Echo command executed, any errors? " + (errCode2 == 0 ? "No" : "Yes"));
            } catch (IOException e) {
                e.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            nrdone.getAndIncrement();

        }
    }


}