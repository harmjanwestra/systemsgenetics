package transeqtlcrossmap;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;


import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.*;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;
import umcg.genetica.util.Primitives;


public class DeterminePolyABias {


//    public static void main(String[] args) {
//        String refgenome = "/Users/harm-jan/SyncData/TMP/polyABias/Homo_sapiens_assembly38.fasta.gz";
//        String gtf = "/Users/harm-jan/SyncData/TMP/polyABias/Homo_sapiens.GRCh38.106.gtf.gz";
//        String outdir = "/Users/harm-jan/SyncData/TMP/polyABias/";
//
//        DeterminePolyABias b = new DeterminePolyABias();
//        try {
//            b.determinePolyABias(refgenome, outdir, gtf);
//        } catch (IOException e) {
//            throw new RuntimeException(e);
//        }
//    }


    protected void determinePolyABias(String genomeFasta, String outputdir, String ensemblannotation) throws IOException {

        GTFAnnotation s = new GTFAnnotation(ensemblannotation);

        // convert gene strings to ensembl gene ids
        ArrayList<Gene> genesToExport = s.getGenesAsArrayList();
        System.out.println(genesToExport.size() + " genes found in annotation.");

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

        TextFile geneStats = new TextFile(outputdir + "allgenes.stats.txt", TextFile.W);

        int blocksize = 5;
        int maxblock = 35;
        String blockStringFull = "";
        String blockStringExonic = "";
        for (int k = blocksize; k < maxblock; k += blocksize) {
            if (k == blocksize) {
                blockStringFull = "NrAStretchesInGeneWithLen" + k;
                blockStringExonic = "NrAStretchesInExonsWithLen" + k;
            } else {
                blockStringFull += "\tNrAStretchesInGeneWithLen" + k;
                blockStringExonic += "\tNrAStretchesInExonsWithLen" + k;
            }
        }

        String header = "Gene\t" +
                "Symbol\t" +
                "Chromosome\t" +
                "Start\t" +
                "Stop\t" +
                "Strand\t" +
                "NrBasesIncludingIntrons\t" +
                "NrA-IncludingIntrons\tNrT-IncludingIntrons\tNrC-IncludingIntrons\tNrG-IncludingIntrons\t" +
                "PropA-IncludingIntrons\tPropT-IncludingIntrons\tPropC-IncludingIntrons\tPropG-IncludingIntrons\t" +
                blockStringFull + "\t" +
                "NrTranscripts\t" +
                "MaxLengthTranscriptName\t" +
                "MaxLengthTranscriptNrExons\t" +
                "MaxLengthTranscriptNrBasesExcludingIntrons\t" +
                "NrA-ExcludingIntrons\tNrT-ExcludingIntrons\tNrC-ExcludingIntrons\tNrG-ExcludingIntrons\t" +
                "PropA-ExcludingIntrons\tPropT-ExcludingIntrons\tPropC-ExcludingIntrons\tPropG-ExcludingIntrons\t" +
                blockStringExonic;
        System.out.println(header);
        geneStats.writeln(header);

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
            int gctr = 0;
            for (Gene g : genesToExport) {
                if (g.getChromosome().equals(seqchr) && seqname.length() < 6) {
                    byte[] genebytes = new byte[g.getSize()];

                    System.arraycopy(bases, g.getStart(), genebytes, 0, g.getSize());
                    String geneStr = byteToStr(genebytes);
//                    String geneATCG = getATCG(geneStr, g.getStrand());

                    // detect polyA stretches

                    int curstretch = 0;

                    Strand strand = g.getStrand();
                    char compchar = 'A';
                    if (strand.equals(Strand.NEG)) {
                        compchar = 'T';
                    }

                    int nrMisMatches = 1;


//                for (Integer key : keysArr) {
//                    System.out.println(key + "\t" + stretchlen.get(key));
//                }
                    // System.out.println();

                    // divide up into blocks of 5
                    HashMap<Integer, Integer> stretchlenFullGene = determinePolyAStretches(geneStr, compchar, nrMisMatches);
                    Set<Integer> keys = stretchlenFullGene.keySet();
                    ArrayList<Integer> keysArr = new ArrayList<>();
                    keysArr.addAll(keys);
                    Collections.sort(keysArr);


                    ArrayList<Integer> blocksfull = new ArrayList<>();
                    for (int k = blocksize; k < maxblock; k += blocksize) {
                        int kctr = 0;
                        for (int ki = 0; ki < keysArr.size(); ki++) {
                            int key = keysArr.get(ki);
                            if (key >= k) {
                                int ct = stretchlenFullGene.get(key);
                                kctr += ct;
                            }
                        }
//                    System.out.println(k + "\t" + kctr);
                        blocksfull.add(kctr);
                    }


                    // System.out.println();
                    // get the exonic regions for the longest transcript
                    Transcript maxTranscript = null;
                    int maxlen = -1;
                    for (Transcript t : g.getTranscripts()) {
                        int len = 0;
                        for (Exon e : t.getExons()) {
                            len += e.getSize();
                        }
//                    System.out.println(t.getName() + "\t" + len);
                        if (len > maxlen) {
                            maxTranscript = t;
                            maxlen = len;
                        }
                    }
//                System.out.println(maxTranscript.getName() + " is max transcript");
//                System.out.println();
                    String geneStringExonic = getBasesAsString(maxTranscript, maxlen, bases);
                    HashMap<Integer, Integer> stretchlenExonic = determinePolyAStretches(geneStringExonic, compchar, nrMisMatches);
                    keys = stretchlenExonic.keySet();
                    keysArr = new ArrayList<>();
                    keysArr.addAll(keys);
                    Collections.sort(keysArr);
                    ArrayList<Integer> blocksexonic = new ArrayList<>();
                    for (int k = blocksize; k < maxblock; k += blocksize) {
                        int kctr = 0;
                        for (int ki = 0; ki < keysArr.size(); ki++) {
                            int key = keysArr.get(ki);
                            if (key >= k) {
                                int ct = stretchlenExonic.get(key);
                                kctr += ct;
                            }
                        }
//                    System.out.println(k + "\t" + kctr);
                        blocksexonic.add(kctr);
                    }


                    String outln = g.getName()
                            + "\t" + g.getGeneSymbol()
                            + "\t" + g.getChromosome()
                            + "\t" + g.getStart()
                            + "\t" + g.getStop()
                            + "\t" + g.getStrand()
                            + "\t" + g.getSize()
                            + "\t" + getATCG(geneStr, strand)
                            + "\t" + Strings.concat(Primitives.toPrimitiveArr(blocksfull), Strings.tab)
                            + "\t" + g.getTranscripts().size()
                            + "\t" + maxTranscript.getName()
                            + "\t" + maxTranscript.getExons().size()
                            + "\t" + maxlen
                            + "\t" + getATCG(geneStringExonic, strand)
                            + "\t" + Strings.concat(Primitives.toPrimitiveArr(blocksexonic), Strings.tab);
                    System.out.println(outln);
                    geneStats.writeln(outln);

//                TextFile gtfo = new TextFile(outputdir + g.getName() + ".txt", TextFile.W);
//                gtfo.write(geneStr);
//                gtfo.close();

                    gctr += 1;
                    if (gctr % 100 == 0) {
                        System.out.println(gctr + " / " + genesToExport.size() + " genes parsed.");
                    }
//                if (gctr == 10) {
//                    geneStats.close();
//                    System.exit(0);
//                }
                }

            }
            System.out.println("Done exporting genes for sequence " + seqname);
            seq = fastaFile.nextSequence();
        }
        geneStats.close();
    }

    private String getATCG(String geneStr, Strand strand) {
        int nrA = 0;
        int nrT = 0;
        int nrC = 0;
        int nrG = 0;
        boolean negStrand = strand.equals(Strand.NEG);
        for (char c : geneStr.toCharArray()) {
            switch (c) {
                case 'A':
                    if (negStrand) {
                        nrT += 1;
                    } else {
                        nrA += 1;
                    }
                    break;
                case 'C':
                    if (negStrand) {
                        nrG += 1;
                    } else {
                        nrC += 1;
                    }
                    break;
                case 'T':
                    if (negStrand) {
                        nrA += 1;
                    } else {
                        nrT += 1;
                    }
                    break;
                case 'G':
                    if (negStrand) {
                        nrC += 1;
                    } else {
                        nrG += 1;
                    }
                    break;
            }
        }

        double total = nrA + nrT + nrC + nrG;
        double propA = 0, propC = 0, propT = 0, propG = 0;
        if (nrA > 0) {
            propA = nrA / total;
        }
        if (nrT > 0) {
            propT = nrT / total;
        }
        if (nrC > 0) {
            propC = nrC / total;
        }
        if (nrG > 0) {
            propG = nrG / total;
        }

        return nrA + "\t" + nrT + "\t" + nrC + "\t" + nrG + "\t" + propA + "\t" + propT + "\t" + propC + "\t" + propG;

    }

    private String getBasesAsString(Transcript transcript, int len, byte[] bases) throws UnsupportedEncodingException {
        byte[] genebytes = new byte[len];
        // System.out.println();
        int curpos = 0;

        // sort exons by position rather than rank
        ArrayList<Exon> tmpExons = new ArrayList<>(transcript.getExons());
        tmpExons.sort(new FeatureComparator(true));

        for (Exon e : tmpExons) {
//            System.out.println(e.getName() + "\t" + e.getStart() + "\t" + transcript.getExonRank(e));
            System.arraycopy(bases, e.getStart(), genebytes, curpos, e.getSize());
            curpos += e.getSize();
        }
        return byteToStr(genebytes);
    }

    private String byteToStr(byte[] bytes) throws UnsupportedEncodingException {
        return new String(bytes, StandardCharsets.UTF_8);
    }

    private HashMap<Integer, Integer> determinePolyAStretches(String geneStr, char compchar, int nrMisMatches) {
        boolean[] isA = new boolean[geneStr.length()];
        char prevbase = 'Q';
        int curstretch = 0;
        HashMap<Integer, Integer> stretchlen = new HashMap<>();

        for (int p = 0; p < geneStr.length(); p++) {
            char curbase = geneStr.charAt(p);
            if (prevbase == compchar) {
                isA[p] = true;
                curstretch += 1;
            } else {

                // check mismatches; is one of the previous bases up to nrMisMatches an A (or a T)?
                boolean aPreviousBaseIsAnA = false;
                if (p > nrMisMatches) {
                    for (int q = p; q > (p - nrMisMatches); q--) {
                        if (isA[q]) {
                            aPreviousBaseIsAnA = true;
                        }
                    }
                }
                if (aPreviousBaseIsAnA) {
                    // add to the current stretch
                    curstretch += 1;
                } else {
                    if (curstretch > 1) {
                        Integer ctr = stretchlen.get(curstretch);
                        if (ctr == null) {
                            ctr = 0;
                        }
                        ctr += 1;
                        stretchlen.put(curstretch, ctr);
                    }
                    curstretch = 0;
                }
            }
            prevbase = curbase;
        }


        return stretchlen;
    }
}