package mbqtl.vcf;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;

import mbqtl.enums.AnalysisType;
import umcg.genetica.features.Feature;
import umcg.genetica.features.SNPFeature;
import umcg.genetica.io.Gpio;


public class VCFVariantProvider implements Iterator<VCFVariant> {

    private final ArrayList<SNPFeature> snpFeaturesForGene;
    private final AnalysisType analysisType;
    private final String origvcfFile;
    private final boolean[] genotypeSamplesToInclude;
    private final Set<String> snpLimitSetForGene;
    private final boolean hasTemplateStr;
    private final String[] vcfFilenames;

    private Iterator<VCFVariant> currentSNPIterator = null;
    private VCFTabix tabix = null;
    private VCFGenotypeData currentVCFGenotypeData = null;
    private int prevChr = -1;
    private int nextStreamingChr = 1;
    private int snpFeatureCounter = 0;
    private boolean debug = false;

    public VCFVariantProvider(Feature region,
                              AnalysisType analysisType,
                              String origvcfFile,
                              boolean[] genotypeSamplesToInclude,
                              Set<String> snpLimitSetForGene,
                              ArrayList<SNPFeature> snpFeaturesForGene) throws IOException {
        this.analysisType = analysisType;
        this.origvcfFile = origvcfFile;
        this.genotypeSamplesToInclude = genotypeSamplesToInclude;
        this.snpLimitSetForGene = snpLimitSetForGene;
        this.snpFeaturesForGene = snpFeaturesForGene;

        this.hasTemplateStr = (origvcfFile.contains("CHR"));
        if (hasTemplateStr) {
            vcfFilenames = new String[23];
            for (int i = 1; i < 23; i++) {
                vcfFilenames[i] = origvcfFile.replace("CHR", "" + i);
                if (!Gpio.exists(vcfFilenames[i])) {
                    vcfFilenames[i] = null;
                }
            }
        } else {
            vcfFilenames = null;
        }

        if (snpFeaturesForGene != null) {
            if (!snpFeaturesForGene.isEmpty()) {
                initTabix(snpFeaturesForGene.get(0));
            }
        } else if (analysisType == AnalysisType.CIS) {
            initTabix(region);
        } else {
            initTrans();
        }


    }

    public void enableDebug() {
        this.debug = true;
    }


//    public Iterator<VCFVariant> provide() {
//
//        Iterator<VCFVariant> snpIterator = null;
//        VCFTabix tabix = null;
//        try {
//            if (analysisType == ANALYSISTYPE.CIS) {
//                // replace vcf file with the one from the actual chromosome, if CHR template has been provided.
//                if (geneChromosomeObj.getNumber() != prevchr && origvcfFile.contains("CHR")) {
//                    // Note: code below is not correct, because it overwrites datasets that are outside of scope of the thread that runs this code
//                    // can cause race conditions!!!!
//                    String actualVCFFile = origvcfFile.replace("CHR", "" + geneChromosomeObj.getNumber());

    /// /						updateDatasets(actualVCFFile); // assume the order of samples is equal across files, for now...
//                    tabix = new VCFTabix(actualVCFFile);
//                    prevchr = geneChromosomeObj.getNumber();
//                } else {
//                    tabix = new VCFTabix(vcfFile);
//                }
//                snpIterator = tabix.getVariants(cisRegion, genotypeSamplesToInclude, snpLimitSetForGene);
//            } else {
//                // TODO: NOT IMPLEMENTED YET
//            }
//        } catch (ArrayIndexOutOfBoundsException e) {
//            System.out.println("Error querying variants for gene: " + gene + " - " + geneSymbol + " - Located on: " + geneChromosomeObj.toString() + " - " + start + " - " + stop);
//        }
//
//    }
    public void close() throws IOException {

        if (tabix != null) {
            tabix.close();
        }
        if (currentVCFGenotypeData != null) {
            currentVCFGenotypeData.close();
        }

    }


    // function for random access
    private void initTabix(Feature queryRegion) throws IOException {

        if (debug) {
            System.out.println("Init tabix, query: " + queryRegion.toString());
        }
        if (!hasTemplateStr) {
            if (debug) {
                System.out.println("VCF init: Option 1 - single file");
            }
//            System.out.println("Option 1");
            if (Gpio.exists(origvcfFile) && tabix == null) {
//                System.out.println("Option 1: open new file handle");
                tabix = new VCFTabix(origvcfFile);
            }
        } else {
            if (debug) {
                System.out.println("VCF init: Option 2 - multiple files");
            }
            int queryChrNumber = queryRegion.getChromosome().getNumber();

            // only open file handle when not already open and variant chr is different from previously queried variant
            if (queryChrNumber != prevChr && queryChrNumber < vcfFilenames.length) {
                String actualVCFFile = vcfFilenames[queryChrNumber];
                if (actualVCFFile != null) {
                    tabix = new VCFTabix(actualVCFFile);
                    prevChr = queryRegion.getChromosome().getNumber();
                } else {
                    tabix = null;
                }
            } else {
                if(debug){
                    System.out.println("VCF init: Option 2 - tabix handle already open.");
                }
            }
//            else {
//                tabix = null;
//            }
            if (debug) {
                System.out.println("Option 2 tabix handle is null for chr " + queryRegion.getChromosome().getNumber() + ": " + (tabix == null) + "\tQuery: " + queryRegion.toString());
            }
        }

        if (tabix != null) {
            if (debug) {
                System.out.println("Tabix: get variants: " + queryRegion);
            }
            currentSNPIterator = tabix.getVariants(queryRegion, genotypeSamplesToInclude, snpLimitSetForGene);
        }
    }

    private void initTrans() throws IOException {
        if (currentVCFGenotypeData != null) {
            currentVCFGenotypeData.close();
            currentVCFGenotypeData = null;
        }

        String vcfFile = origvcfFile;
        if (hasTemplateStr) {
            if (nextStreamingChr < vcfFilenames.length) {
                vcfFile = vcfFilenames[nextStreamingChr];
                nextStreamingChr++;
            } else {
                vcfFile = null;
            }
        }

        if (vcfFile != null && Gpio.exists(vcfFile)) {
            currentVCFGenotypeData = new VCFGenotypeData(vcfFile, genotypeSamplesToInclude, snpLimitSetForGene);
        }
    }

    @Override
    public boolean hasNext() {
        try {
            if (snpFeaturesForGene != null) {
//                System.out.println(snpFeatureCounter);
                if (debug) {
                    System.out.println("HasNext, current ctr:" + snpFeatureCounter);
                }
                return snpFeatureCounter < snpFeaturesForGene.size();
            } else if (analysisType == AnalysisType.CIS) {
                if (currentSNPIterator == null) {
                    return false;
                }
                return currentSNPIterator.hasNext();
            } else {

                if (currentVCFGenotypeData == null) {
                    return false;
                }

                boolean currentDataHasNext = currentVCFGenotypeData.hasNext();
                if (currentDataHasNext) {
                    return true;
                } else if (hasTemplateStr) {
                    initTrans();
                    return currentVCFGenotypeData.hasNext();
                }
                return false;
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public VCFVariant next() {
        if (this.hasNext()) {
            if (snpFeaturesForGene != null) {

                // iterate variant list for gene until we've found a match by random lookup in the VCF using tabix
                while (snpFeatureCounter < snpFeaturesForGene.size()) {
                    try {
                        SNPFeature feature = snpFeaturesForGene.get(snpFeatureCounter);
                        if (debug) {
                            System.out.println("Looking for feature: " + snpFeatureCounter + "\t" + feature.getName());
                        }
                        initTabix(feature);
                        if (debug) {
                            System.out.println("Current snp iterator is not null: " + ((currentSNPIterator != null)));
                        }
                        if (currentSNPIterator != null) {
                            if (debug) {
                                System.out.println("Current snp iterator has next: " + (currentSNPIterator.hasNext()));
                            }
                            while (currentSNPIterator.hasNext()) {
                                VCFVariant variant = currentSNPIterator.next();
                                if (debug) {
                                    if (variant != null) {
                                        System.out.println("VCF has variant: " + snpFeatureCounter + "\t" + variant.getId() + "\t" + feature.toString());
                                    } else {
                                        System.out.println("VCF returns null variant for: " + snpFeatureCounter + "\t" + feature.getName() + "\t" + feature.toString());
                                    }
                                }
//
                                if (variant != null && variant.getId().equals(feature.getName())) {
                                    if (debug) {
                                        System.out.println("Matching ID: " + variant);
                                        System.out.println("Feature ctr: " + snpFeatureCounter);
                                        System.out.println(variant.getHwep());
                                        System.out.println(variant.getMAF());
                                        System.out.println(variant.getCallrate());
                                    }
                                    snpFeatureCounter++;
                                    return variant;
                                }
                            }
                        }
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }

                    snpFeatureCounter++;
                    if (debug) {
                        System.out.println("Adding to snpfeaturectr: " + snpFeatureCounter);
                    }
                }
            } else if (analysisType == AnalysisType.CIS) {
                return currentSNPIterator.next();

            } else {
                return currentVCFGenotypeData.next();
            }
        }

        return null;
    }
}
