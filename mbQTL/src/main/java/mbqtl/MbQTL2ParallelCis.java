package mbqtl;

import mbqtl.datastructures.Dataset;
import mbqtl.stat.BetaDistributionMLE;
import mbqtl.stat.FisherWeightedMetaAnalysis;
import mbqtl.stat.PVal;
import mbqtl.stat.RankArray;
import mbqtl.vcf.VCFTabix;
import mbqtl.vcf.VCFVariant;
import umcg.genetica.console.ProgressBar;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.enums.Strand;
import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.Descriptives;
import umcg.genetica.math.stats.Heterogeneity;
import umcg.genetica.math.stats.ZScores;
import umcg.genetica.text.Strings;
import umontreal.iro.lecuyer.probdist.BetaDist;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class MbQTL2ParallelCis extends QTLAnalysis {


    DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
    private DecimalFormat dfDefault = new DecimalFormat("#.######", symbols);
    private DecimalFormat dfPval = new DecimalFormat("#.####E0", symbols);
    private int nrPermutations = 1000;
    private int cisWindow = 1000000;
    private long randomSeed = 123456789;
    private boolean rankData = true;
    private boolean outputAll = false;
    private boolean outputSNPLog = false;
    private boolean replaceMissingGenotypes = false;
    private boolean dumpPermutationPvalues = false;
    private ANALYSISTYPE analysisType = ANALYSISTYPE.CIS;
    private boolean testNonParseableChr = false;
    private METAANALYSISMETHOD metaanalysismethod = METAANALYSISMETHOD.EMP;

    public MbQTL2ParallelCis(String vcfFile, int chromosome, String linkfile, String snpLimitFile, String geneLimitFile, String snpGeneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, String outfile) throws IOException {
        super(vcfFile, chromosome, linkfile, snpLimitFile, geneLimitFile, snpGeneLimitFile, geneExpressionDataFile, geneAnnotationFile, outfile);
        if (datasets.length < minNumberOfDatasets) {
            System.out.println(minNumberOfDatasets + " datasets required, but only " + datasets.length + " datasets defined. Changing setting to: " + datasets.length);
            minNumberOfDatasets = datasets.length;
        }
    }

    public void setNrPermutations(int nrPermutations) {
        this.nrPermutations = nrPermutations;
    }

    public void setCisWindow(int cisWindow) {
        this.cisWindow = cisWindow;
    }

    public void setRandomSeed(long randomSeed) {
        this.randomSeed = randomSeed;
    }

    public void setRankData(boolean rankData) {
        this.rankData = rankData;
    }

    public void setOutputAll(boolean outputAll) {
        this.outputAll = outputAll;
    }

    public void setOutputAllPermutations(boolean dumpPermutationPvalues) {
        this.dumpPermutationPvalues = dumpPermutationPvalues;
    }

    public void setOutputSNPLog(boolean outputSNPLog) {
        this.outputSNPLog = outputSNPLog;
    }

    public void setReplaceMissingGenotypes(boolean replaceMissingGenotypes) {
        this.replaceMissingGenotypes = replaceMissingGenotypes;
    }

    public void setMetaanalysismethod(METAANALYSISMETHOD method) {
        this.metaanalysismethod = method;
    }

    public void run() throws IOException {
        /*
        TODO:
        x function to test specific SNPs
        x output dataset Z-scores and sample sizes
        x plotting?
        x check whether code runs with no permutations
        - a way to run on a VCF containing all chromosomes
        - check proper way to shuffle when there is missing data in the genotype as well as the phenotype data
        */

        System.out.println();
        System.out.println("----------------------------------------------------------");
        System.out.println("QTL meta-analysis with Beta distribution approximated null");
        System.out.println("----------------------------------------------------------");
        System.out.println("MAF:\t" + mafthreshold);
        System.out.println("Call-rate:\t" + callratethreshold);
        System.out.println("HWE-P:\t" + hwepthreshold);
        System.out.println("Min nr datasets:\t" + minNumberOfDatasets);
        System.out.println("Nr Permutations:\t" + nrPermutations);
        System.out.println("Cis window:\t" + cisWindow);
        System.out.println("Random randomSeed:\t" + randomSeed);
        System.out.println("Ranking data:\t" + rankData);
//        System.out.println("Re-ranking data:\t");
        System.out.println("Replacing missing genotypes: " + replaceMissingGenotypes);
        System.out.println("Min observations: " + minObservations);
        System.out.println("Min genotype count:\t" + minGenotypeCount);
        System.out.println("Splitting multi allelic variants:\t" + splitMultiAllelics);
        System.out.println("Skipping multi allelic variants:\t" + (!splitMultiAllelics));
        System.out.println("Writing all snp/feature pairs: " + outputAll);
        System.out.println("Writing SNP log: " + outputSNPLog);
        System.out.println("Writing all permutations: " + dumpPermutationPvalues);
        System.out.println("Meta-analysis method:\t" + metaanalysismethod);
        System.out.println();

        Chromosome chromosomeObj = Chromosome.parseChr("" + chromosome);
        System.out.println("Processing: " + vcfFile);
//		if (!Gpio.exists(vcfFile)) {
//
//		}

        // initialize permutation seeds
        long[] randomSeed = new long[nrPermutations];
        Random rand = new Random(this.randomSeed);
        for (int i = 0; i < randomSeed.length; i++) {
            randomSeed[i] = rand.nextLong();
        }

        String datasetstr = "";
        for (int d = 0; d < datasets.length; d++) {
            if (d == 0) {
                datasetstr = datasets[d].getName();
            } else {
                datasetstr += ";" + datasets[d].getName();
            }
        }

        // initialize output

        TextFile outTopFx = new TextFile(outputPrefix + "-TopEffects.txt", TextFile.W);
        String[] headerTopFx = new String[26];
        headerTopFx[0] = "Gene";
        headerTopFx[1] = "GeneChr";
        headerTopFx[2] = "GenePos";
        headerTopFx[3] = "GeneStrand";
        headerTopFx[4] = "GeneSymbol";
        headerTopFx[5] = "SNP";
        headerTopFx[6] = "SNPChr";
        headerTopFx[7] = "SNPPos";
        headerTopFx[8] = "SNPAlleles";
        headerTopFx[9] = "SNPEffectAllele";
        headerTopFx[10] = "SNPEffectAlleleFreq";
        headerTopFx[11] = "MetaP";
        headerTopFx[12] = "MetaPN";
        headerTopFx[13] = "MetaPZ";
        headerTopFx[14] = "MetaBeta";
        headerTopFx[15] = "MetaSE";
        headerTopFx[16] = "MetaI2";
        headerTopFx[17] = "NrDatasets";
        headerTopFx[18] = "DatasetCorrelationCoefficients(" + datasetstr + ")";
        headerTopFx[19] = "DatasetZScores(" + datasetstr + ")";
        headerTopFx[20] = "DatasetSampleSizes(" + datasetstr + ")";
        headerTopFx[21] = "NrTestedSNPs";
        headerTopFx[22] = "ProportionBetterPermPvals";
        headerTopFx[23] = "BetaDistAlpha";
        headerTopFx[24] = "BetaDistBeta";
        headerTopFx[25] = "BetaAdjustedMetaP";
        if (metaanalysismethod == METAANALYSISMETHOD.FISHERZFIXED || metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
            headerTopFx[14] = "MetaR";
            headerTopFx[19] = "DatasetFisherZ(" + datasetstr + ")";
            if (metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                headerTopFx[11] += "-Random";
                headerTopFx[14] += "-Random";
                headerTopFx[19] += "-Random";
            }
        }
        String headerTopFxStr = Strings.concat(headerTopFx, Strings.tab);
        if (geneGroups != null) {
            headerTopFxStr = "Group\t" + headerTopFxStr;
        }
        outTopFx.writeln(headerTopFxStr);
        TextFile outAll = null;
        if (outputAll) {
            outAll = new TextFile(outputPrefix + "-AllEffects.txt.gz", TextFile.W);
//            String headerAll = "Gene\tGeneSymbol\tSNP\tSNPAlleles\tSNPEffectAllele\tMetaP\tMetaPN\tMetaPZ\tMetaBeta\tMetaSE\tNrDatasets\tProportionBetterPermPvals\tBetaAdjustedMetaP";
            String[] headerAll = new String[21];
            headerAll[0] = "Gene";
            headerAll[1] = "GeneChr";
            headerAll[2] = "GenePos";
            headerAll[3] = "GeneStrand";
            headerAll[4] = "GeneSymbol";
            headerAll[5] = "SNP";
            headerAll[6] = "SNPChr";
            headerAll[7] = "SNPPos";
            headerAll[8] = "SNPAlleles";
            headerAll[9] = "SNPEffectAllele";
            headerAll[10] = "SNPEffectAlleleFreq";
            headerAll[11] = "MetaP";
            headerAll[12] = "MetaPN";
            headerAll[13] = "MetaPZ";
            headerAll[14] = "MetaBeta";
            headerAll[15] = "MetaSE";
            headerAll[16] = "MetaI2";
            headerAll[17] = "NrDatasets";
            headerAll[18] = "DatasetCorrelationCoefficients(" + datasetstr + ")";
            headerAll[19] = "DatasetZScores(" + datasetstr + ")";
            headerAll[20] = "DatasetSampleSizes(" + datasetstr + ")";
            if (metaanalysismethod == METAANALYSISMETHOD.FISHERZFIXED || metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                headerAll[14] = "MetaR";
                headerAll[19] = "DatasetFisherZ(" + datasetstr + ")";
                if (metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                    headerAll[11] += "-Random";
                    headerAll[14] += "-Random";
                    headerAll[19] += "-Random";
                }
            }
            String headerAllStr = Strings.concat(headerAll, Strings.tab);
            if (geneGroups != null) {
                headerAllStr = "Group\t" + headerAllStr;
            }
            outAll.writeln(headerAllStr);
        }

        TextFile permutationoutput = null;
        if (dumpPermutationPvalues) {
            permutationoutput = new TextFile(outputPrefix + "-Permutations.txt.gz", TextFile.W);
            String permHeader = "";
            for (int i = 0; i < nrPermutations; i++) {
                if (i == 0) {
                    permHeader += "Perm" + i;
                } else {
                    permHeader += "\tPerm" + i;
                }
            }
            if (geneGroups != null) {
                permutationoutput.writeln("Group\tGene\tSNP\t" + permHeader);
            } else {
                permutationoutput.writeln("Gene\tSNP\t" + permHeader);
            }
        }

        // iterate genes
        TextFile finalOutAll = outAll;
        AtomicInteger testedgenes = new AtomicInteger();
        TextFile logout = new TextFile(outputPrefix + "-log.txt.gz", TextFile.W);
        TextFile snplogout = null;
        if (outputSNPLog) {
            snplogout = new TextFile(outputPrefix + "-snpqclog.txt.gz", TextFile.W);
            String snplogheader = "Gene\tVariant\tRefAllele/AltAllele\tOVerallAltAlleleFreq\tNrTotalAlleles\tNrDatasetsPassingQC\tPassQC(" + datasetstr + ")\tMAF(" + datasetstr + ")\tCR(" + datasetstr + ")\tHWE-P(" + datasetstr + ")";
            if (geneGroups != null) {
                snplogheader = "Group\t" + snplogheader;
            }
            snplogout.writeln(snplogheader);
        }

        TextFile finalSnplogout = snplogout;
        TextFile finalOutAll1 = outAll;
        TextFile finalPermutationoutput = permutationoutput;
        if (geneGroups != null) {
            ProgressBar pb = new ProgressBar(geneGroups.size(), "Processing " + geneGroups.size() + " groups of genes...");

            geneGroups.entrySet().parallelStream().forEach(entry -> {
                String groupName = entry.getKey();
                HashSet<String> groupGenes = entry.getValue();
                try {
                    testGenes(groupName, groupGenes, logout, finalSnplogout, finalOutAll1, finalPermutationoutput, outTopFx, randomSeed, testedgenes);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                pb.iterateSynchedPrint();
            });
            pb.close();
        } else {
            ArrayList<Gene> geneObjs = new ArrayList<>();
            for (int g = 0; g < expressionData.genes.length; g++) {
                Gene obj = geneAnnotation.getGene(expressionData.genes[g]);
                if (obj != null) {
                    geneObjs.add(obj);
                }
            }
            Collections.sort(geneObjs, new FeatureComparator());
            ProgressBar pb = new ProgressBar(expressionData.genes.length, "Processing " + expressionData.genes.length + " genes...");
            IntStream.range(0, geneObjs.size()).parallel().forEach(g -> {
                HashSet<String> groupGenes = new HashSet<>();
                groupGenes.add(geneObjs.get(g).getName());
                try {
                    testGenes(null, groupGenes, logout, finalSnplogout, finalOutAll1, finalPermutationoutput, outTopFx, randomSeed, testedgenes);
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
//			pb.set(g + 1);
                pb.iterateSynchedPrint();
            }); // ENDIF: iterate genes
            pb.close();
        }

        outTopFx.close();
        if (outputAll) {
            outAll.close();
        }
        if (dumpPermutationPvalues) {
            permutationoutput.close();
        }

        logout.close();
        if (outputSNPLog) {
            snplogout.close();
        }
        TextFile outFinished = new TextFile(outputPrefix + "-TopEffects.finished", TextFile.W);
        outFinished.writeln("Tested genes:\t" + testedgenes.get());
        outFinished.close();
    }

    private void testGenes(String groupName,
                           HashSet<String> genes,
                           TextFile logout,
                           TextFile snpLogOut,
                           TextFile outAll,
                           TextFile outPermutation,
                           TextFile outTopFx,
                           long[] seed,
                           AtomicInteger testedgenes) throws IOException {

        double[] permutationPvals = new double[nrPermutations];
        Arrays.fill(permutationPvals, 1);
        final UnpermutedResult topUnpermutedResult = new UnpermutedResult(); // this is safe, because values are only changed once per SNP, when permutation == -1

        int nrTestsPerformed = 0;

        // todo: check if the genes in the group are in the same genomic location
        // if so, variants can be cached, which will save some disk and parsing overhead
        HashSet<String> uniqueSNPIDs = new HashSet<>();
        int prevchr = -1;
        for (String gene : genes) {
//            Integer geneAnnotationId = geneAnnotation.getGeneId(gene);
            Gene geneAnnotationObj = geneAnnotation.getGene(gene);
            if (geneAnnotationObj == null) {
                logout.writelnsynced("Skipping " + gene + " since it has no annotation.");
            } else {
                int g = expressionData.geneMap.get(gene);
                double[] expData = expressionData.data[g];

                // define CIS window
                int pos = geneAnnotationObj.getStart(); // .getStartPos(geneAnnotationId);
                Strand strand = geneAnnotationObj.getStrand(); // .getStrand(geneAnnotationId);
                String geneSymbol = geneAnnotationObj.getGeneSymbol(); // .getSymbol(geneAnnotationId);
                int start = pos - cisWindow;
                if (start < 0) {
                    start = 0;
                }
                int stop = pos + cisWindow;
//				int chr = geneAnnotationObj.getChromosome().getNumber(); // .getChr(geneAnnotationId);
                Chromosome geneChromosomeObj = geneAnnotationObj.getChromosome();

                if (geneChromosomeObj.getNumber() < 1 && !testNonParseableChr) {
                    logout.writelnsynced("Skipping " + gene + " since maps to an unparseable chromosome.");
                    continue;
                }

                Feature cisRegion = new Feature(geneChromosomeObj, start, stop);

                // split expression data per dataset
                double[][] expressionPerDataset = new double[datasets.length][];
                IntStream.range(0, datasets.length).forEach(d -> {
                    Dataset thisDataset = datasets[d];
                    double[] datasetExp = thisDataset.select(expData, thisDataset.getExpressionIds());
                    double[] datasetExpRanked = datasetExp;

                    if (rankData) {
                        RankArray ranker = new RankArray();
                        datasetExpRanked = ranker.rank(datasetExp, true); // does this work with NaNs? answer: no
                    }
//                        for (int v = 0; v < datasetExp.length; v++) {
//                            System.out.println(thisDataset.name + "\t" + datasetExp[v] + "\t" + datasetExpRanked[v]);
//                        }
                    expressionPerDataset[d] = datasetExpRanked;
                });
//                    System.exit(0);

                Set<String> snpLimitSetForGene = snpLimitSet;
                if (snpGeneLimitSet != null) {
                    snpLimitSetForGene = snpGeneLimitSet.get(gene);
                }


                Iterator<VCFVariant> snpIterator = null;
                VCFTabix tabix = null;
                try {
                    if (analysisType == ANALYSISTYPE.CIS) {
                        // replace vcf file with the one from the actual chromosome, if CHR template has been provided.
                        if (geneChromosomeObj.getNumber() != prevchr && origvcfFile.contains("CHR")) {
                            // Note: code below is not correct, because it overwrites datasets that are outside of scope of the thread that runs this code
                            // can cause race conditions!!!!
                            String actualVCFFile = origvcfFile.replace("CHR", "" + geneChromosomeObj.getNumber());
//						updateDatasets(actualVCFFile); // assume the order of samples is equal across files, for now...
                            tabix = new VCFTabix(actualVCFFile);
                            prevchr = geneChromosomeObj.getNumber();
                        } else {
                            tabix = new VCFTabix(vcfFile);
                        }
                        snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSetForGene);
                    } else {
                        // TODO: NOT IMPLEMENTED YET
                    }
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.out.println("Error querying variants for gene: " + gene + " - " + geneSymbol + " - Located on: " + geneChromosomeObj.toString() + " - " + start + " - " + stop);
                }

                // TODO: can the VCF parsing be cached somehow?

                while (snpIterator != null && snpIterator.hasNext()) {
                    VCFVariant variantTmp = snpIterator.next();


                    ArrayList<VCFVariant> variants;
                    if (variantTmp != null) {
//                        String variantId2 = variantTmp.getId();
//                        if (variantId2.equals("rs779116082;rs745987535")) {
//                            System.out.println("Got it!");
//                            splitMultiAllelics = true;
//                        }

                        if (!variantTmp.isMultiallelic()) {
                            variants = new ArrayList<>();
                            variants.add(variantTmp);
                        } else if (splitMultiAllelics) {
                            variants = variantTmp.splitMultiAllelic();
                        } else {
                            variants = new ArrayList<>();
                        }
                    } else {
                        variants = new ArrayList<>();
                    }

                    for (VCFVariant variant : variants) {
                        String variantId = variant.getId();
                        uniqueSNPIDs.add(variantId);
                        if ((snpLimitSet == null || snpLimitSet.contains(variantId)) ||
                                (snpGeneLimitSet == null || (snpGeneLimitSet.containsKey(gene) && snpGeneLimitSet.get(gene).contains(variantId)))
                        ) {
                            final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
                            double[] dosages = null;
                            if (useHardGenotypeCalls) {
                                dosages = getGenotype(variant.getGenotypesAsByteVector());
                            } else {
                                dosages = getDosage(variant.getDosage());
                            }


                            // split genotype data per dataset, perform QC
                            double[][] genotypesPerDataset = new double[datasets.length][];
                            double[][] dosagesPerDataset = new double[datasets.length][];
                            VariantQCObj[] qcobjs = new VariantQCObj[datasets.length];
                            double[] finalDosages = dosages;
                            IntStream.range(0, datasets.length).forEach(d -> {
                                Dataset thisDataset = datasets[d];
//                                if (dosages == null || genotypes == null) {
//                                    System.out.println("Error? " + variant.getId() + "\t" + variants.size());
//                                }
                                dosagesPerDataset[d] = thisDataset.select(finalDosages, thisDataset.getGenotypeIds()); // select required dosages
                                genotypesPerDataset[d] = thisDataset.select(genotypes, thisDataset.getGenotypeIds()); // select required genotype IDs

                                VariantQCObj qcobj = checkVariant(genotypesPerDataset[d]);
                                if (qcobj.passqc) {
                                    if (replaceMissingGenotypes) {
                                        // only replace missing genotypes on variants that pass the qc thresholds
                                        double meanDosage = Util.meanGenotype(dosagesPerDataset[d]);
                                        double meanGenotype = Util.meanGenotype(genotypesPerDataset[d]);
                                        for (int i = 0; i < dosagesPerDataset[d].length; i++) {
                                            if (genotypesPerDataset[d][i] == -1) {
                                                genotypesPerDataset[d][i] = meanGenotype;
                                                dosagesPerDataset[d][i] = meanDosage;
                                            }
                                        }
                                    }


                                    // prune the data here once, to check if there are enough values to go ahead with this snp/gene combo
                                    // but only if the variant is passing the QC in the first place for the samples selected in this dataset
                                    Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(genotypesPerDataset[d],
                                            dosagesPerDataset[d],
                                            expressionPerDataset[d]);

                                    // check the variant again, taking into account missingness in the expression data
                                    qcobj = checkVariant(prunedDatasetData.getLeft());

                                    // require minimum number of observations, otherwise kick out dataset from analysis
                                    if (prunedDatasetData.getLeft().length < minObservations) {
                                        qcobj.passqc = false;
                                    }
                                }
                                qcobjs[d] = qcobj;
                            });

//                                System.out.println("");
                            // run permutations, and non-permuted result (permutation == -1)
                            AtomicBoolean tested = new AtomicBoolean(false);
//							TextFile finalOutAll = outAll;
                            double[] permutationPvalsForSNP = null;
                            if (dumpPermutationPvalues) {
                                permutationPvalsForSNP = new double[nrPermutations];
                                Arrays.fill(permutationPvalsForSNP, 1);
                            }
                            double[] finalPermutationPvalsForSNP = permutationPvalsForSNP;
                            IntStream.range(-1, nrPermutations).forEach(permutation -> {
                                FisherWeightedMetaAnalysis fisherZ = null;

                                if (metaanalysismethod == METAANALYSISMETHOD.FISHERZFIXED || metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                                    fisherZ = new FisherWeightedMetaAnalysis();
                                }
                                double[] zscores = new double[datasets.length];
                                double[] correlations = new double[datasets.length];
                                int[] samplesizes = new int[datasets.length];
                                Arrays.fill(samplesizes, -1);
                                Arrays.fill(zscores, Double.NaN);
                                Arrays.fill(correlations, Double.NaN);

                                // iterate datasets
                                int nrAltAlleles = 0;
                                int nrTotalAlleles = 0;
                                int dsWithMinObs = 0;
                                int nrsnpspassqc = 0;

                                for (int d = 0; d < datasets.length; d++) {
                                    Dataset thisDataset = datasets[d];
                                    double[] datasetGt = genotypesPerDataset[d]; // thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs
                                    VariantQCObj qcobj = qcobjs[d]; // check maf, hwep, call-rate, number of genotypes per genotype group
                                    if (qcobj.passqc) {
                                        nrsnpspassqc++;
                                        double[] datasetExp = expressionPerDataset[d];
                                        double[] datasetExpCopy = new double[datasetExp.length];
                                        System.arraycopy(datasetExp, 0, datasetExpCopy, 0, datasetExpCopy.length);

                                        double[] datasetDs = dosagesPerDataset[d];

                                        // if this is a permutation, shuffle the data
                                        if (permutation != -1) {
                                            Util.shuffleArray(datasetExpCopy, seed[permutation]);
                                        }

                                        // prune the data (remove missing values)
                                        // can't prune the data earlier (would save a lot of compute time) because shuffling is performed over all available samples for this dataset
                                        // this is because the order of permuted samples should be equal across all SNPs
                                        Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetGt,
                                                datasetDs,
                                                datasetExpCopy);

                                        // re-rank data here? original EMP does not, but it is the right thing to do...
                                        double[] datasetExpPruned = prunedDatasetData.getRight();
                                        if (datasetExpPruned.length >= minObservations) {

//                                    if (rankData) {
//                                        RankArray ranker = new RankArray();
//                                        datasetExpPruned = ranker.rank(datasetExpPruned, true); // does this work with NaNs? answer: no
//                                    }


                                            double[] datasetDsPruned = prunedDatasetData.getMiddle();
//											double[] datasetDsPrunedCopy = prunedDatasetData.getMiddle();
                                            double[] datasetGtPruned = prunedDatasetData.getLeft();
//											datasetDsPruned = prunedDatasetData.getMiddle();

                                            // datasetGtPruned = Util.centerScale(prunedDatasetData.getLeft());

//                                                for (int v = 0; v < datasetDsPruned.length; v++) {
//                                                    System.out.println(thisDataset.name + "\t" + datasetDsPruned[v] + "\t" + datasetExpPruned[v]);
//                                                }

                                            // perform correlation
                                            if (Descriptives.variance(datasetDsPruned) > 0 && Descriptives.variance(datasetExpPruned) > 0) {

                                                dsWithMinObs++;
                                                // count the number of alleles, used later to estimate Beta and SE from MetaZ
                                                if (permutation == -1) {
                                                    for (int i = 0; i < datasetGtPruned.length; i++) {

                                                        if (datasetDsPruned[i] >= 0.5 && datasetDsPruned[i] <= 1.5) {
                                                            nrAltAlleles += 1;
                                                        } else if (datasetDsPruned[i] > 1.5) {
                                                            nrAltAlleles += 2;
                                                        }
//														System.out.println(d + "\t" + i + "\t" + datasetGtPruned[i] + "\t" + datasetDsPruned[i] + "\t" + nrAltAlleles);
                                                    }
                                                    nrTotalAlleles += datasetGtPruned.length * 2;
                                                }

                                                datasetDsPruned = Util.centerScale(datasetDsPruned); // This may not be required
                                                datasetExpPruned = Util.centerScale(datasetExpPruned);

                                                double r = Correlation.correlate(datasetDsPruned, datasetExpPruned);
//												double r2 = Correlation.correlate(datasetDsPruned, datasetExpPruned, 0d, 0d, 1d, 1d);
//												double deltar = Math.abs(r2) - Math.abs(r);
//												System.out.println(r + "\t" + r2 + "\t" + deltar);
                                                double p = PVal.getPvalue(r, datasetExpPruned.length - 2);
                                                double z = ZScores.pToZTwoTailed(p); // p value is already two-tailed, so need to use this other p-value conversion method... :/; returns negative z-scores by default
                                                if (r > 0) {
                                                    z *= -1; // flip z-score if correlation is positive because p-value conversion returns only negative z-scores
                                                }

                                                // prevent edge-cases
                                                if (!Double.isNaN(r)) {
                                                    zscores[d] = z;
                                                    correlations[d] = r;
                                                    samplesizes[d] = datasetExpPruned.length;
                                                }
                                            }
                                        } // endif nrobservations >= minobservations
                                    } // endif qcobj.passqc
                                } // ENDfor: test every dataset

                                // determine number of datasets with data
                                int nDatasets = 0;
                                int totalSampleSize = 0;
                                for (int d = 0; d < zscores.length; d++) {
                                    if (!Double.isNaN(zscores[d])) {
                                        totalSampleSize += samplesizes[d];
                                        nDatasets++;
                                    }
                                }
//								System.out.println();
//								System.out.println(variantId + "\t" + nrAltAlleles + "\t" + nrTotalAlleles);
                                double overallAltAlleleFreq = (double) nrAltAlleles / nrTotalAlleles;
                                // write some log stuff
                                if (permutation == -1) {
                                    try {
                                        if (outputSNPLog) {
                                            String passStr = "";
                                            String mafStr = "";
                                            String crStr = "";
                                            String hweStr = "";
                                            for (int d = 0; d < qcobjs.length; d++) {
                                                if (d == 0) {
                                                    if (qcobjs[d].passqc) {
                                                        passStr = "T";
                                                    } else {
                                                        passStr = "F";
                                                    }
                                                    mafStr = "" + dfDefault.format(qcobjs[d].maf);
                                                    crStr = "" + dfDefault.format(qcobjs[d].cr);
                                                    hweStr = "" + toNeatP(qcobjs[d].hwep);
                                                } else {
                                                    if (qcobjs[d].passqc) {
                                                        passStr += ";T";
                                                    } else {
                                                        passStr += ";F";
                                                    }
                                                    mafStr += ";" + dfDefault.format(qcobjs[d].maf);
                                                    crStr += ";" + dfDefault.format(qcobjs[d].cr);
                                                    hweStr += ";" + toNeatP(qcobjs[d].hwep);
                                                }
                                            }

                                            String snplogStr = gene + "\t" +
                                                    variantId + "\t" +
                                                    variant.getAlleles()[0] + "/" +
                                                    variant.getAlleles()[1] + "\t" +
                                                    dfDefault.format(overallAltAlleleFreq) + "\t" +
                                                    nrTotalAlleles + "\t" +
                                                    nrsnpspassqc + "\t" +
                                                    passStr + "\t" +
                                                    mafStr + "\t" +
                                                    crStr + "\t" +
                                                    hweStr;
                                            if (geneGroups != null) {
                                                snplogStr = groupName + "\t" + snplogStr;
                                            }
                                            snpLogOut.writelnsynced(snplogStr);
                                        }

                                        if (snpGeneLimitSet != null && snpGeneLimitSet.get(gene) != null) {
                                            logout.writelnsynced(gene + "\t" + variantId + " effect is present in  " + nDatasets + " and has " + dsWithMinObs + " with >= " + minObservations + ", " + nrsnpspassqc + " of the dataset SNPs pass QC thresholds");
                                        }
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                }

                                if (nDatasets >= minNumberOfDatasets) {
                                    // meta-analyze, weight by sample size
                                    double metaZ = 0;
                                    double metaP = 1;
                                    double metaBeta = 0;
                                    double metaBetaSE = 1;
                                    double metaI2 = 100;

                                    double[] zScoresToPrint = zscores;
                                    if (metaanalysismethod == METAANALYSISMETHOD.FISHERZFIXED || metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                                        double[] metaZtmp;
                                        if (metaanalysismethod == METAANALYSISMETHOD.FISHERZRANDOM) {
                                            metaZtmp = fisherZ.metaAnalyzeCorrelationsRandomEffect(correlations, samplesizes);
                                        } else {
                                            metaZtmp = fisherZ.metaAnalyzeCorrelationsFixedEffect(correlations, samplesizes);
                                        }
                                        metaZ = metaZtmp[0];
                                        metaP = metaZtmp[1];
                                        metaBeta = metaZtmp[2];
                                        metaBetaSE = metaZtmp[3];
                                        metaI2 = metaZtmp[6];
                                        for (int d = 0; d < correlations.length; d++) {
                                            double r = correlations[d];
                                            double z = zscores[d];
                                            if (!Double.isNaN(z)) {
                                                zScoresToPrint[d] = fisherZ.correlationToFisherZ(r);
                                            }
                                        }
                                    } else {
                                        metaZ = ZScores.getWeightedZ(zscores, samplesizes);
                                        metaP = ZScores.zToP(metaZ);
                                        if (permutation == -1) {
                                            double[] betaAndSEEstimate = ZScores.zToBeta(metaZ, overallAltAlleleFreq, totalSampleSize);
                                            metaBeta = betaAndSEEstimate[0];
                                            metaBetaSE = betaAndSEEstimate[1];
                                            Triple<Double, Double, Integer> metaIsq = Heterogeneity.getISq(zscores, samplesizes);
                                            metaI2 = metaIsq.getLeft();
                                        }
                                    }

                                    // calculate heterogeneity
                                    if (permutation != -1) { // this is a permuted result
                                        if (metaP < permutationPvals[permutation]) {
                                            permutationPvals[permutation] = metaP;
                                        }
                                        if (dumpPermutationPvalues) {
                                            finalPermutationPvalsForSNP[permutation] = metaP;
                                        }
                                    } else { // this is a non-permuted result
                                        tested.getAndSet(true);

                                        // calculate overall MAF

//                                        double[] betaAndSEEstimate = ZScores.zToBeta(metaZ, overallAltAlleleFreq, totalSampleSize);
//                                        Triple<Double, Double, Integer> metaIsq = Heterogeneity.getISq(zscores, samplesizes);
//                                        double metaI2 = metaIsq.getLeft();
                                        // non-permuted p-value
                                        if (metaP <= topUnpermutedResult.metaP) {
                                            boolean replace = true;
                                            // if the SNP is in perfect LD (has equal pvalue), select the closest one to the gene
                                            if (metaP == topUnpermutedResult.metaP && topUnpermutedResult.snpID != null) {
                                                // if the Z-score is sufficiently large (>40) we exceed the range of the normal distribution, returning a p-value of ~2x10-232
                                                // in that case, compare the absolute Z-scores to determine the top effect for this gene
                                                if (Math.abs(metaZ) < Math.abs(topUnpermutedResult.metaPZ)) {
                                                    replace = false;
                                                } else if (Math.abs(metaZ) > Math.abs(topUnpermutedResult.metaPZ)) {
                                                    replace = true;
                                                } else { // if the Z-scores are also equal (unlikely)
                                                    int genePos = geneAnnotationObj.getStart(); // .getStartPos(geneAnnotationId);
                                                    int tssDist = Math.abs(genePos - variant.getPos());
                                                    int tssDist2 = Math.abs(genePos - topUnpermutedResult.snpPos);
                                                    if (tssDist > tssDist2) {
                                                        replace = false;
                                                    }
                                                }
                                            }
                                            if (replace) {
                                                topUnpermutedResult.gene = gene;
                                                topUnpermutedResult.metaP = metaP;
                                                topUnpermutedResult.metaPN = totalSampleSize;
                                                topUnpermutedResult.metaPZ = metaZ;
                                                topUnpermutedResult.metaPD = nDatasets;
                                                topUnpermutedResult.metaI2 = metaI2;
                                                topUnpermutedResult.zscores = zScoresToPrint;
                                                topUnpermutedResult.samplesizes = samplesizes;
                                                topUnpermutedResult.correlations = correlations;
                                                topUnpermutedResult.snpEffectAlleleFreq = overallAltAlleleFreq;
                                                topUnpermutedResult.metaBeta = metaBeta;
                                                topUnpermutedResult.metaBetaSE = metaBetaSE;
                                                topUnpermutedResult.snpID = variant.getId();
                                                topUnpermutedResult.snpChr = variant.getChr();
                                                topUnpermutedResult.snpPos = variant.getPos();
                                                topUnpermutedResult.snpAlleles = variant.getAlleles()[0] + "/" + variant.getAlleles()[1];
                                                topUnpermutedResult.snpEffectAllele = variant.getAlleles()[1];
                                            }
                                        }

                                        if (outputAll) { // this code only runs when in the 'not-permuted' iteration
                                            String snpAlleles = variant.getAlleles()[0] + "/" + variant.getAlleles()[1];
                                            String snpEffectAllele = variant.getAlleles()[1];

                                            String outln = gene
                                                    + "\t" + geneChromosomeObj.getNumber()
                                                    + "\t" + pos
                                                    + "\t" + strand
                                                    + "\t" + geneSymbol
                                                    + "\t" + variant.getId()
                                                    + "\t" + variant.getChr()
                                                    + "\t" + variant.getPos()
                                                    + "\t" + snpAlleles
                                                    + "\t" + snpEffectAllele
                                                    + "\t" + dfDefault.format(overallAltAlleleFreq)
                                                    + "\t" + metaP
                                                    + "\t" + totalSampleSize
                                                    + "\t" + dfDefault.format(metaZ)
                                                    + "\t" + dfDefault.format(metaBeta)
                                                    + "\t" + dfDefault.format(metaBetaSE)
                                                    + "\t" + dfDefault.format(metaI2)
                                                    + "\t" + nDatasets
                                                    + "\t" + toNeatStr(correlations)
                                                    + "\t" + toNeatStr(zScoresToPrint)
                                                    + "\t" + toNeatSampleSizeStr(samplesizes);
                                            if (geneGroups != null) {
                                                outln = groupName + "\t" + outln;
                                            }
                                            try {
                                                outAll.writelnsynced(outln);
                                            } catch (IOException e) {
                                                e.printStackTrace();
                                            }

                                        }
                                    }
                                }
                            }); // ENDIF: for each permutation


                            // determine if SNP was tested somehow...
                            if (tested.get()) {
                                if (dumpPermutationPvalues) {
                                    try {
                                        String outln = gene + "\t" + variant.getId() + "\t" + Strings.concat(finalPermutationPvalsForSNP, Strings.tab);
                                        if (geneGroups != null) {
                                            outln = groupName + "\t" + outln;
                                        }
                                        outPermutation.writelnsynced(outln);
                                    } catch (IOException e) {
                                        e.printStackTrace();
                                    }
                                }
                                nrTestsPerformed++;
                            }
                        }  // ENDIF: variant in snpGeneLimitSet || snpLimitSet
                    } // ENDIF: if variant != null
                } // ENDIF: while snpiterator has next
                if (tabix != null) {
                    tabix.close();
                }
            } // end if annotation not null
        } // end iterate group's genes

        if (nrTestsPerformed == 0) {
            if (geneGroups != null) {
                logout.writelnsynced(groupName + " has " + uniqueSNPIDs.size() + " SNPs in the CIS-window, but none passed QC.");
            } else {
                String gene = genes.toArray(new String[0])[0];
                logout.writelnsynced(gene + " has " + uniqueSNPIDs.size() + " SNPs in the CIS-window, but none passed QC.");
            }
        } else {
            // determine beta distribution etc
            double propBetterPvals = 0;
            double betaAdjPval = 1;
            double[] shape = new double[]{Double.NaN, Double.NaN};
            BetaDist betaDistribution = null;
            boolean output = true;
            if (nrPermutations > 1) { // permutations are required for the following step
                for (int p = 0; p < permutationPvals.length; p++) {
                    if (permutationPvals[p] <= topUnpermutedResult.metaP) {
                        propBetterPvals++;
                    }
                }
                propBetterPvals /= permutationPvals.length;  // permutation p-value (equals 0.0 if metaP is very small)
                try {
                    BetaDistributionMLE mle = new BetaDistributionMLE();
                    shape = mle.fit(permutationPvals);
                    betaDistribution = new BetaDist(shape[0], shape[1]);
                    betaAdjPval = betaDistribution.cdf(topUnpermutedResult.metaP);
                    if (betaAdjPval < 2.0E-323D) {
                        betaAdjPval = 2.0E-323D;
                    }
                } catch (org.apache.commons.math3.exception.TooManyEvaluationsException tmee) {
                    if (geneGroups != null) {
                        logout.writelnsynced(groupName + " failed: Beta MLE Model did not converge.");
                    } else {
                        String gene = genes.toArray(new String[0])[0];
                        logout.writelnsynced(gene + " failed: Beta MLE Model did not converge.");
                    }

                    betaAdjPval = propBetterPvals;
                    output = false;
                }
            } else {
                propBetterPvals = 1;
            }

            if (output) {
                testedgenes.getAndIncrement();
//				Integer geneAnnotationId = geneAnnotation.getGeneId(topUnpermutedResult.gene);
                Gene geneAnnotationObj = geneAnnotation.getGene(topUnpermutedResult.gene);
                int pos = geneAnnotationObj.getStart(); // .getStartPos(geneAnnotationId);
                Strand strand = geneAnnotationObj.getStrand(); // .getStrand(geneAnnotationId);
                String geneSymbol = geneAnnotationObj.getGeneSymbol(); // .getSymbol(geneAnnotationId);
                String outln = topUnpermutedResult.gene
                        + "\t" + geneAnnotationObj.getChromosome().getNumber()
                        + "\t" + pos + "\t" + strand + "\t" + geneSymbol
                        + "\t" + topUnpermutedResult.snpID
                        + "\t" + topUnpermutedResult.snpChr
                        + "\t" + topUnpermutedResult.snpPos
                        + "\t" + topUnpermutedResult.snpAlleles
                        + "\t" + topUnpermutedResult.snpEffectAllele
                        + "\t" + dfDefault.format(topUnpermutedResult.snpEffectAlleleFreq)
                        + "\t" + topUnpermutedResult.metaP
                        + "\t" + topUnpermutedResult.metaPN
                        + "\t" + dfDefault.format(topUnpermutedResult.metaPZ)
                        + "\t" + dfDefault.format(topUnpermutedResult.metaBeta)
                        + "\t" + dfDefault.format(topUnpermutedResult.metaBetaSE)
                        + "\t" + dfDefault.format(topUnpermutedResult.metaI2)
                        + "\t" + topUnpermutedResult.metaPD
                        + "\t" + toNeatStr(topUnpermutedResult.correlations)
                        + "\t" + toNeatStr(topUnpermutedResult.zscores)
                        + "\t" + toNeatSampleSizeStr(topUnpermutedResult.samplesizes)
                        + "\t" + nrTestsPerformed // not sure if this is correct for groups with > 1 genes in them
                        + "\t" + dfDefault.format(propBetterPvals)
                        + "\t" + dfDefault.format(shape[0])
                        + "\t" + dfDefault.format(shape[1])
                        + "\t" + betaAdjPval;
                if (geneGroups != null) {
                    outln = groupName + "\t" + outln;
                }
                outTopFx.writelnsynced(outln);
            }
        }
    }

    private String toNeatP(double pval) {
        if (pval <= 0) {
            return "0";
        } else if (pval == 1) {
            return "1";
        } else if (pval < 0.00001) {
            return "" + dfPval.format(pval);
        } else {
            return "" + dfDefault.format(pval);
        }
    }

    private String toNeatStr(double[] vals) {
        StringBuilder b = new StringBuilder();
        for (int d = 0; d < vals.length; d++) {
            if (d > 0) {
                b.append(";");
            }
            if (Double.isNaN(vals[d])) {
                b.append("-");
            } else {
                b.append("" + dfDefault.format(vals[d]));
            }
        }
        return b.toString();
    }

    private String toNeatStr(int[] vals) {
        StringBuilder b = new StringBuilder();
        for (int d = 0; d < vals.length; d++) {
            if (d > 0) {
                b.append(";");
            }
            b.append("" + vals[d]);
        }
        return b.toString();
    }

    private String toNeatSampleSizeStr(int[] vals) {
        StringBuilder b = new StringBuilder();
        for (int d = 0; d < vals.length; d++) {
            if (d > 0) {
                b.append(";");
            }
            if (vals[d] < 0) {
                b.append("-");
            } else {
                b.append("" + vals[d]);
            }
        }
        return b.toString();
    }

    public void setTestNonParseableChr() {
        this.testNonParseableChr = true;
    }

    enum METAANALYSISMETHOD {
        EMP,
        FISHERZFIXED,
        FISHERZRANDOM,
    }

    enum ANALYSISTYPE {
        CIS,
        TRANS,
        CISTRANS
    }

    private class UnpermutedResult {
        public double metaBeta;
        public double metaBetaSE;
        public String snpID;
        public String snpAlleles;
        public String snpEffectAllele;
        public String snpChr;
        public int snpPos;
        public double snpEffectAlleleFreq;
        public double[] zscores;
        public int[] samplesizes;
        public double[] correlations;
        public String gene;
        public double metaI2 = 1;
        double metaP = 1;
        double metaPN = 0;
        double metaPZ = 0;
        double metaPD = 0;

        @Override
        public String toString() {
            return "UnpermutedResult{" +
                    "metaP=" + metaP +
                    ", metaPN=" + metaPN +
                    ", metaPZ=" + metaPZ +
                    ", metaPD=" + metaPD +
                    ", snpID=" + snpID +
                    ", gene=" + gene +
                    '}';
        }

    }


}
