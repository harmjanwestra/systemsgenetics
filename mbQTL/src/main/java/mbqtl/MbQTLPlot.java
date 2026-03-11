package mbqtl;

import JSci.maths.ArrayMath;
import com.itextpdf.text.DocumentException;
import mbqtl.datastructures.Dataset;
import mbqtl.enums.MetaAnalysisMethod;
import mbqtl.gfx.ForestPlotEQTL;
import mbqtl.gfx.ForestplotPanel;
import mbqtl.gfx.QTLPanel;
import mbqtl.stat.FisherWeightedMetaAnalysis;
import mbqtl.stat.PVal;
import mbqtl.stat.RankArray;
import mbqtl.vcf.VCFVariant;
import mbqtl.vcf.VCFVariantProvider;
import umcg.genetica.containers.Triple;
import umcg.genetica.enums.Chromosome;
import umcg.genetica.features.Feature;
import umcg.genetica.features.FeatureComparator;
import umcg.genetica.features.Gene;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.graphics.Grid;
import umcg.genetica.math.stats.Correlation;
import umcg.genetica.math.stats.ZScores;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import java.util.stream.IntStream;


public class MbQTLPlot extends QTLAnalysis {

    private int nrPermutations = 1000;
    private int cisWindow = 1000000;
    private long randomSeed = 123456789;
    private boolean rankData = true;
    private boolean outputAll = false;
    private boolean replaceMissingGenotypes = false;
    private MetaAnalysisMethod metaanalysismethod = MetaAnalysisMethod.EMP;
    private boolean dontrankvariants;

    public MbQTLPlot(String vcfFile, int chromosome, String linkfile, String snpLimitFile, String geneLimitFile, String snpGeneLimitFile, String geneExpressionDataFile, String geneAnnotationFile, int minNumberOfDatasets,
                     int minObservations, String outputPrefix) throws IOException {
        super(vcfFile, chromosome, linkfile, snpLimitFile, geneLimitFile, snpGeneLimitFile, geneExpressionDataFile, geneAnnotationFile, minNumberOfDatasets, minObservations, outputPrefix);
    }

    public void setNrPermutations(int nrPermutations) {
        this.nrPermutations = nrPermutations;
    }

    public void setMetaanalysismethod(MetaAnalysisMethod method) {
        this.metaanalysismethod = method;
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

    public void setReplaceMissingGenotypes(boolean replaceMissingGenotypes) {
        this.replaceMissingGenotypes = replaceMissingGenotypes;
    }

    public void plot() throws IOException {
        System.out.println();
        System.out.println("-----------");
        System.out.println("QTL plotter");
        System.out.println("-----------");


        for (int g = 0; g < expressionData.genes.length; g++) {
            String gene = expressionData.genes[g];
            Gene geneAnnotationObj = geneAnnotation.getGene(gene);

            if (geneAnnotationObj != null) {
                // get chromosome
                int genechr = geneAnnotationObj.getChromosome().getNumber();
                Chromosome chromosomeObj = geneAnnotationObj.getChromosome();
                System.out.println(g + "/" + expressionData.genes.length + " - Plotting " + gene + " on chr " + genechr);


//                String tmpVcf = origvcfFile;
//                if (origvcfFile.contains("CHR")) {
//                    tmpVcf = origvcfFile.replaceAll("CHR", "" + chr);
//                }
//                System.out.println("Opening " + tmpVcf);
//                VCFTabix tabix = new VCFTabix(tmpVcf);


                double[] expData = expressionData.data[g];

                // define CIS window
                int genepos = geneAnnotationObj.getStart(); // .getStartPos(geneAnnotationId);
                String geneSymbol = geneAnnotationObj.getGeneSymbol(); //.getSymbol(geneAnnotationId);
                int cisstart = genepos - cisWindow;
                if (cisstart < 0) {
                    cisstart = 0;
                }
                int cisstop = genepos + cisWindow;
                Feature cisRegion = new Feature(chromosomeObj, cisstart, cisstop);
                System.out.println("Region: " + cisRegion.toString());

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
                    expressionPerDataset[d] = datasetExpRanked;

                });

                System.out.println("Querying tabix..");
                Set<String> snpLimitSetForGene = snpLimitSet;
                ArrayList<SNPFeature> snpLimitSetForGeneFeatures = null;
                if (snpGeneLimitSet != null) {
                    snpLimitSetForGene = snpGeneLimitSet.get(gene);
                }
                if (snpLimitSetForGene != null && snpAnnotation != null) {
                    snpLimitSetForGeneFeatures = new ArrayList<>();
                    for (String snp : snpLimitSetForGene) {
                        Integer snpid = snpAnnotation.getId(snp);
                        if (snpid != null) {
                            int snppos = snpAnnotation.getPos(snpid);
                            Chromosome snpchr = snpAnnotation.getChr(snpid);
                            SNPFeature f = new SNPFeature(snpchr, snppos, snppos + 1);
                            f.setName(snp);
                            snpLimitSetForGeneFeatures.add(f);
                        }
                    }
                    snpLimitSetForGeneFeatures.sort(new FeatureComparator()); // sort by position
//                    System.out.println();
//                    System.out.println(snpLimitSetForGeneFeatures.size() + " SNPs for gene: " + gene);
                }
                VCFVariantProvider snpIterator = new VCFVariantProvider(cisRegion, analysisType, origvcfFile, genotypeSamplesToInclude, snpLimitSetForGene, snpLimitSetForGeneFeatures);

//                Iterator<VCFVariant> snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSet);

                while (snpIterator.hasNext()) {
                    VCFVariant variant = snpIterator.next();

                    if (variant != null) {
                        String variantId = variant.getId();
                        if ((snpLimitSet == null || snpLimitSet.contains(variantId)) || (snpGeneLimitSet == null || (snpGeneLimitSet.containsKey(gene) && snpGeneLimitSet.get(gene).contains(variantId)))) {
                            final double[] genotypes = getGenotype(variant.getGenotypesAsByteVector());
                            final double[] dosages = getDosage(variant.getDosage());

                            // split genotype data per dataset, perform QC
                            double[][] genotypesPerDataset = new double[datasets.length][];
                            double[][] dosagesPerDataset = new double[datasets.length][];
                            VariantQCObj[] qcobjs = new VariantQCObj[datasets.length];
                            IntStream.range(0, datasets.length).forEach(d -> {
                                Dataset thisDataset = datasets[d];
                                dosagesPerDataset[d] = thisDataset.select(dosages, thisDataset.getGenotypeIds()); // select required dosages
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
                                    Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(genotypesPerDataset[d], dosagesPerDataset[d], expressionPerDataset[d]);

                                    // check the variant again, taking into account missingness in the expression data
                                    qcobj = checkVariant(prunedDatasetData.getLeft());

                                    // require minimum number of observations, otherwise kick out dataset from analysis
                                    if (prunedDatasetData.getLeft().length < minObservations) {
                                        qcobj.passqc = false;
                                    }
                                }
                                qcobjs[d] = qcobj;
                            });

                            int nrPassingQC = 0;
                            for (int d = 0; d < qcobjs.length; d++) {
                                if (qcobjs[d].passqc) {
                                    nrPassingQC++;
                                }
                            }
                            if (nrPassingQC >= minNumberOfDatasets) {
//                                Grid grid = new Grid(200, 200, 1, nrPassingQC, 100, 100);
                                Grid qtlPlotGrid = new Grid(300, 300, datasets.length/4, datasets.length / 4, 100, 50);
//                                QTLPanel metaAnalysisPanel = new QTLPanel(1, 1);
                                ForestplotPanel metaAnalysisPanel = new ForestplotPanel(2,2);
                                if(metaanalysismethod == MetaAnalysisMethod.EMP){
                                    metaAnalysisPanel.setStatistic(ForestplotPanel.STATISTIC.ZSCORE);
                                } else if(metaanalysismethod == MetaAnalysisMethod.FISHERZFIXED || metaanalysismethod == MetaAnalysisMethod.FISHERZRANDOM){
                                    metaAnalysisPanel.setStatistic(ForestplotPanel.STATISTIC.CORRELATION);
                                }

                                ArrayList<ForestPlotEQTL> qtls = new ArrayList<>();
                                ArrayList<QTLPanel> panels = new ArrayList<>();

                                // iterate datasets
                                String geneStr = gene + " (" + geneSymbol + ")";
//                                ArrayList<Double> metaValsX = new ArrayList<>();
//                                ArrayList<Double> metaValsY = new ArrayList<>();

                                double[] zscores = new double[nrPassingQC];

                                double[] correlations = new double[nrPassingQC];
                                int[] samplesizes = new int[nrPassingQC];

                                int passqcctr = 0;
                                int sumN = 0;
                                for (int d = 0; d < datasets.length; d++) {
                                    Dataset thisDataset = datasets[d];
                                    double[] datasetGt = genotypesPerDataset[d]; // thisDataset.select(genotypes, thisDataset.genotypeIds); // select required genotype IDs
                                    VariantQCObj qcobj = qcobjs[d]; // check maf, hwep, call-rate, number of genotypes per genotype group
                                    if (!qcobj.passqc) {
                                        QTLPanel dsPanel = new QTLPanel(1, 1);
                                        dsPanel.setNotTested();
                                        dsPanel.setDatasetDetails(thisDataset.getName(), geneStr, variantId, 0, 1, 0, 0, 0);
                                        panels.add(dsPanel);
                                    } else {
                                        double[] datasetExp = expressionPerDataset[d];
                                        double[] datasetExpCopy = new double[datasetExp.length];
                                        System.arraycopy(datasetExp, 0, datasetExpCopy, 0, datasetExpCopy.length);

                                        double[] datasetDs = dosagesPerDataset[d];
//                                        System.out.println(gene);
//                                        for (int q = 0; q < datasetGt.length; q++) {
//                                            System.out.println(thisDataset.name + "\t" + datasetGt[q] + "\t" + datasetDs[q] + "\t" + datasetExp[q]);
//                                        }

                                        // prune the data (remove missing values)
                                        // can't prune the data earlier (would save a lot of compute time) because shuffling is performed over all available samples for this dataset
                                        // this is because the order of permuted samples should be equal across all SNPs
                                        Triple<double[], double[], double[]> prunedDatasetData = pruneMissingValues(datasetGt, datasetDs, datasetExpCopy);

                                        // re-rank data here? original EMP does not, but it is the right thing to do...
                                        double[] datasetExpPruned = prunedDatasetData.getRight();
//                                    if (rankData) {
//                                        RankArray ranker = new RankArray();
//                                        datasetExpPruned = ranker.rank(datasetExpPruned, true); // does this work with NaNs? answer: no
//                                    }
                                        double[] origExp = datasetExpPruned;
                                        double[] datasetDsPruned = prunedDatasetData.getMiddle();

//                                        datasetExpPruned = Util.centerScale(datasetExpPruned);
                                        double[] datasetGtPruned = prunedDatasetData.getLeft();
                                        if(useHardGenotypeCalls){
                                            datasetDsPruned = datasetGtPruned;
                                        }

                                        boolean passgtcount = true;
                                        if(minGenotypeCount>0){
                                            int aa = 0;
                                            int ab = 0;
                                            int bb = 0;
                                            for(int v=0;v<datasetDsPruned.length;v++){
                                                if(datasetDsPruned[v] < 0.5){
                                                    aa++;
                                                } else if(datasetDsPruned[v]>1.5){
                                                    bb++;
                                                } else {
                                                    ab++;
                                                }
                                            }

                                            if(aa < minGenotypeCount || ab < minGenotypeCount || bb < minGenotypeCount){
                                                passgtcount = false;
                                            }
                                        }

                                        if (rankData) {
                                            RankArray ranker = new RankArray();
                                            datasetExpPruned = ranker.rank(datasetExpPruned, true);

                                            // also rank genotypes if no hard calls are used
                                            if(!useHardGenotypeCalls){
                                                if(!dontrankvariants) {
                                                    datasetDsPruned = ranker.rank(datasetDsPruned, true);
                                                }
                                            }
                                        }


                                        if (datasetExpPruned.length == 0) {
                                            System.out.println(datasets[d].getName() + ": exp: " + datasetExpCopy.length + "\tgt: " + datasetDs.length + "\t has no genotypes after pruning, but passqc: " + qcobj.passqc);
                                            System.out.println(datasets[d].getName() + ": cr: " + qcobj.cr + "\thwep: " + qcobj.hwep + "\tmaf: " + qcobj.maf);
                                            QTLPanel dsPanel = new QTLPanel(1, 1);
                                            dsPanel.setNotTested();
                                            dsPanel.setDatasetDetails(thisDataset.getName(), geneStr, variantId, 0, 1, 0, 0, 0);
                                            panels.add(dsPanel);
                                        } else {

                                            QTLPanel dsPanel = new QTLPanel(1, 1);
                                            if(!passgtcount){
                                                dsPanel.setNotTested();
                                                zscores[passqcctr] = Double.NaN;
                                                correlations[passqcctr] = Double.NaN;
                                                samplesizes[passqcctr] = 0;
                                                dsPanel.setDatasetDetails(thisDataset.getName(), geneStr, variantId, 0, 1, 0, 0, 0);
                                            } else {
                                                // perform correlation
                                                double r = Correlation.correlate(datasetDsPruned, datasetExpPruned);
                                                double p = PVal.getPvalue(r, datasetExpPruned.length - 2);
                                                double z = ZScores.pToZTwoTailed(p); // p value is already two-tailed, so need to use this other p-value conversion method... :/; returns negative z-scores by default
                                                if (r > 0) {
                                                    z *= -1; // flip z-score if correlation is positive because p-value conversion returns only negative z-scores
                                                }

                                                zscores[passqcctr] = z;
                                                correlations[passqcctr] = r;
                                                samplesizes[passqcctr] = datasetExpPruned.length;
                                                sumN += datasetDsPruned.length;
                                                System.out.println(r+"\t"+p);

                                                if(rankData){
                                                    if(!useHardGenotypeCalls){
                                                        if(!dontrankvariants) {
                                                            // scale back the ranked genotypes to between 0 and 2
                                                            // ranked data starts at 0
                                                            double max = ArrayMath.max(datasetDsPruned);
                                                            double min = ArrayMath.min(datasetDsPruned);
                                                            double delta = max - min;

                                                            System.out.println("Min: " + min + "\tmax: " + max);
                                                            for (int v = 0; v < datasetDsPruned.length; v++) {
                                                                datasetDsPruned[v] = ((datasetDsPruned[v] - min) / delta) * 2;
                                                            }
                                                        }
                                                    }
                                                }

                                                ForestPlotEQTL qtl = new ForestPlotEQTL();
                                                qtl.EA = variant.getAlleles()[1];
                                                qtl.snp = variantId;
                                                qtl.dataset = thisDataset.getName();
                                                qtl.samplesize = datasetExpPruned.length;
                                                qtl.isMeta=false;

                                                if (metaanalysismethod == MetaAnalysisMethod.FISHERZFIXED || metaanalysismethod == MetaAnalysisMethod.FISHERZRANDOM) {
                                                    qtl.effectsize = r;
                                                    double rsq = r * r;
                                                    int n = datasetDsPruned.length;
                                                    if(n < 3){
                                                        qtl.error = 0;
                                                    } else {
                                                        qtl.error = (1 - rsq) / (Math.sqrt(n - 3)); // Bonett 2008 approximation of standard error
                                                    }
                                                    qtl.effecttype = ForestPlotEQTL.EFFECTTYPE.CORRELATION;
                                                } else {
                                                    qtl.effectsize = z;
                                                    qtl.error = Double.MIN_VALUE; // no error for Zscore
                                                    qtl.effecttype = ForestPlotEQTL.EFFECTTYPE.ZSCORE;
                                                }

                                                qtl.pvalue = p;
                                                qtls.add(qtl);

                                                dsPanel.setData(datasetDsPruned, datasetExpPruned);
                                                dsPanel.setAlleles(variant.getAlleles());
                                                int n = datasetDsPruned.length;
                                                double maf = qcobj.maf;
                                                dsPanel.setDatasetDetails(thisDataset.getName(), geneStr, variantId, z, p, r, maf, n);
                                            }

                                            panels.add(dsPanel);
                                            passqcctr++;
                                        }


                                    }
                                }

                                // define meta-analysis panel
                                if(datasets.length>1) {
//                                    metaAnalysisPanel.setData(Primitives.toPrimitiveArr(metaValsX), Primitives.toPrimitiveArr(metaValsY));
//                                    metaAnalysisPanel.setAlleles(variant.getAlleles());

                                    ForestPlotEQTL qtl = new ForestPlotEQTL();
                                    if (metaanalysismethod == MetaAnalysisMethod.FISHERZFIXED || metaanalysismethod == MetaAnalysisMethod.FISHERZRANDOM) {
                                        FisherWeightedMetaAnalysis fisherZ = new FisherWeightedMetaAnalysis();
                                        double[] metaZtmp;
                                        if (metaanalysismethod == MetaAnalysisMethod.FISHERZRANDOM) {
                                            metaZtmp = fisherZ.metaAnalyzeCorrelationsRandomEffect(correlations, samplesizes);
                                            qtl.dataset = "Meta-analysis (Fisher-random)";
                                        } else {
                                            metaZtmp = fisherZ.metaAnalyzeCorrelationsFixedEffect(correlations, samplesizes);
                                            qtl.dataset = "Meta-analysis  (Fisher-fixed)";
                                        }

                                        double metaZ = metaZtmp[0];
                                        double metaP = metaZtmp[1];
                                        double metaR = metaZtmp[2];
                                        double metaBetaSE = metaZtmp[3];
                                        double metaI2 = metaZtmp[6];
                                        qtl.effectsize = metaR;
                                        qtl.effecttype = ForestPlotEQTL.EFFECTTYPE.CORRELATION;
                                        qtl.error =metaBetaSE;
                                        qtl.pvalue = metaP;
                                    } else {
                                        double metaZ = ZScores.getWeightedZ(zscores, samplesizes);
                                        double metaP = ZScores.zToP(metaZ);
                                        qtl.effectsize = metaZ;
                                        qtl.error = Double.MIN_VALUE;
                                        qtl.pvalue = metaP;
                                        qtl.dataset = "Meta-Analysis (Samplesize weighted)";
                                        qtl.effecttype = ForestPlotEQTL.EFFECTTYPE.ZSCORE;
                                    }

                                    qtl.isMeta = true;
                                    qtl.snp = variantId;

                                    qtl.EA = variant.getAlleles()[1];
                                    qtl.samplesize = sumN;
                                    qtls.add(qtl);
                                    metaAnalysisPanel.setData(qtls);

                                    qtlPlotGrid.addPanel(metaAnalysisPanel);

                                }

                                for (QTLPanel p : panels) {
                                    qtlPlotGrid.addPanel(p);
                                }

                                gene = gene.replaceAll("\\|", "_");
                                gene = gene.replaceAll("\\.", "-");
                                geneSymbol = geneSymbol.replaceAll("\\|", "_");
                                geneSymbol = geneSymbol.replaceAll("\\.", "-");
                                if (gene.length() > 20) {
                                    gene = gene.substring(0, 20);
                                }
                                String snp = variant.getId().replaceAll(":", "_");
                                try {
                                    //                                    String fileout = outputPrefix + "-" + gene + "-" + snp + ".pdf";
                                    //                                    System.out.println("Plotting: " + fileout);
                                    //                                    grid.draw(fileout);
                                    String fileout = outputPrefix + "qtlplot-" + genechr + "_" + geneAnnotationObj.getStart() + "-" + gene + "-" + geneSymbol + "_" + snp + ".pdf";
                                    System.out.println("Plotting: " + fileout);
                                    qtlPlotGrid.draw(fileout);

                                } catch (DocumentException e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    public void setDontRankGenotypes() {
        this.dontrankvariants = true;

    }
}
