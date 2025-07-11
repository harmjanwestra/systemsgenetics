package mbqtl;

import mbqtl.enums.AnalysisType;
import mbqtl.enums.MetaAnalysisMethod;
import mbqtl.enums.PermutationStrategy;
import org.apache.commons.cli.*;

import java.io.IOException;

public class Main {

    public static void main(String[] args) {

        Options options = new Options();

        options.addRequiredOption("m", "mode", true, "Mode: [metaqtl|mbqtl|mbqtlsingleds|mbqtlplot|regressqtl|sortfile|determineld|determineldgwas|concatconditional]");
        options.addOption("v", "vcf", true, "Tabix indexed VCF");
        options.addOption("e", "exp", true, "Expression matrix (can be gzipped)");
        options.addOption("eg", "expgroups", true, "File defining groups of phenotypes");
        options.addOption("chr", "chr", true, "Chromosome number");
        options.addOption("g", "gte", true, "Genotype to expression to dataset linkfile (tab separated)");
        options.addOption("sgl", "snpgenelimit", true, "SNP-gene limit file (one line per snp-gene combination, tab separated)");
        options.addOption("sl", "snplimit", true, "SNP limit file (one line per gene ID)");
        options.addOption("gl", "genelimit", true, "Gene limit file (one line per gene ID)");
        options.addOption("a", "annotation", true, "Gene annotation file");
        options.addOption("o", "out", true, "Output prefix");
        options.addOption("seed", "seed", true, "Random seed [default: 123456789]");
        options.addOption("perm", "perm", true, "Number of permutations [default: 1000]");
        options.addOption("ciswindow", "ciswindow", true, "Cis window size [default: 1mb]");
        options.addOption("cis", "cis", false, "Perform cis-QTL analysis (default)");
        options.addOption("trans", "trans", false, "Perform trans-QTL analysis (exclude variants in cis-window)");
        options.addOption("cistrans", "cistrans", false, "Perform cis+trans-QTL analysis (test all variants)");
        options.addOption("maf", "maf", true, "Minor allele frequency threshold [default: 0.01]");
        options.addOption("cr", "cr", true, "Call-rate threshold [default: 0.95]");
        options.addOption("hwep", "hwep", true, "Hardy-Weinberg p-value threshold [default: 0.0001]");
        options.addOption("snpannotation", "snpannotation", true, "SNP annotation file, tab separated: SNPID chr pos");
        options.addOption("minobservations", "minobservations", true, "Require at least this many observations per dataset (i.e. non-NaN genotypes/phenotypes) [default: 10]");
        options.addOption("norank", "norank", false, "Do not rank expression data");
        options.addOption("outputall", "outputall", false, "Output all associations, not just top association per gene");
        options.addOption("outputallpermutations", "outputallpermutations", false, "Output all permuted associations, not just top association per gene");
        options.addOption("snplog", "snplog", false, "Output SNP summary stats per snp/gene pair.");
        options.addOption("nrdatasets", "nrdatasets", true, "Minimum number of datasets required in meta-analysis [default: 2]");
        options.addOption("replacemissinggenotypes", "replacemissinggenotypes", false, "Replace missing genotypes with average genotype: use this when both genotypes and expression data have missing values and perm > 0");
        options.addOption("input", "input", true, "Input file");
        options.addOption("input2", "input2", true, "Input file 2");
        options.addOption("nriters", "nriters", true, "Nr iters to concatenate using concatconditional");
        options.addOption("sortbyz", "sortbyz", false, "Sort by Z-score");
        options.addOption("sortbysnppos", "sortbysnppos", false, "Sort by SNP position");
        options.addOption("eqtlset", "eqtlset", true, "[determineldgwas] - List of eQTL snps to test - txt.gz file");
        options.addOption("eqtlfile", "eqtlfile", true, "[determineldgwas] - eQTL file - txt.gz file");
        options.addOption("gwasset", "gwasset", true, "[determineldgwas] - List of GWAS snps to test - txt.gz file");
        options.addOption("gwasassocfile", "gwasassocfile", true, "[determineldgwas] - GWAS association file - txt.gz file");
        options.addOption("gwaslistfile", "gwaslistfile", true, "[determineldgwas] - GWAS study file - txt.gz file");
        options.addOption("prunedistance", "prunedistance", true, "[determineldgwas] - Pruning distance [default: 1mb]");
        options.addOption("prunethreshold", "prunethreshold", true, "[determineldgwas] - Pruning LD threshold [default: 0.2]");
        options.addOption("ldthreshold", "ldthreshold", true, "[determineldgwas] - LD threshold [default: 0.8]");
        options.addOption("skipchr6", "skipchr6", false, "[determineldgwas] - Skip chr6 when calculating LD");
        options.addOption("matchbyrsid", "matchbyrsid", false, "[determineldgwas] - Match by RsId in stead of full variant id");
        options.addOption("testnonparseablechr", "testnonparseablechr", false, "[mbqtl] - Test variants and genes that map to non-oarseable chromosomes (e.g. patch chromosomes)");
        options.addOption("empzmeta", "empzmeta", false, "Perform fixed effects meta-analysis using legacy eQTL mapping pipeline meta-analysis method [default]");
        options.addOption("fisherzmeta", "fisherzmeta", false, "Perform fixed effects meta-analysis using weighted fisher Z method");
        options.addOption("fisherzmetarandom", "fisherzmetarandom", false, "Perform random effects meta-analysis using weighted fisher Z method");
        options.addOption("mingenotypecount", "mingenotypecount", true, "Minimal number of individuals per genotype group [default: 0]");
        options.addOption("splitmultiallelic", "splitmultiallelic", false, "Split multi allelic variants [default: skip multi allelic variants].");
        options.addOption("usehardgenotypecalls", "usehardgenotypecalls", false, "Use hard genotype calls in stead of dosages derived from genotype probabilities. [default: use dosages if available]");
        options.addOption("onlytestsnps", "onlytestsnps", false, "Only test SNPs, skip indels, STRs, etc. [default: test all variants in VCF passing QC]");
        options.addOption("correlationweights", "correlationweights", true, "Weight the correlation for these weights. Format: sample, weight, tab-separated");
        options.addOption("permutationstrategy", "permutationstrategy", true, "Permutation strategy [matched|random|semirandom]: use same permuted order for each variant and gene per permutation (matched), use completely random assignments (random), or somewhere in between (semirandom). Only random and semirandom can be used if there is missingness in both genotype and genotype data.");

        try {
            CommandLineParser parser = new DefaultParser();
            CommandLine cmd = parser.parse(options, args);

            String mode = cmd.getOptionValue("mode");


            int minNrDatasets = 2;
            if (cmd.hasOption("nrdatasets")) {
                minNrDatasets = Integer.parseInt(cmd.getOptionValue("nrdatasets"));

            }

            int minObservations = 10;
            if (cmd.hasOption("minobservations")) {
                minObservations = Integer.parseInt(cmd.getOptionValue("minobservations"));
            }

            String vcf = null;
            if (cmd.hasOption("vcf")) {
                vcf = cmd.getOptionValue("vcf");
            }

            int chrom = -1;
            if (cmd.hasOption("chr")) {
                chrom = Integer.parseInt(cmd.getOptionValue("chr"));
            }

            String linkfile = null;
            if (cmd.hasOption("gte")) {
                linkfile = cmd.getOptionValue("gte");
            }

            String genelimit = null;
            if (cmd.hasOption("genelimit")) {
                genelimit = cmd.getOptionValue("genelimit");
            }

            String snplimit = null;
            if (cmd.hasOption("snplimit")) {
                snplimit = cmd.getOptionValue("snplimit");
            }
            String snpgenelimit = null;
            if (cmd.hasOption("snpgenelimit")) {
                if (!cmd.hasOption("snpannotation")) {
                    System.out.println("Suggestion: you can also specify --snpannotation when using --snpgenelimit for additional performance if you are interested in few snps/variants.");
                }
                snpgenelimit = cmd.getOptionValue("snpgenelimit");
            }

            String genexpression = null;
            if (cmd.hasOption("exp")) {
                genexpression = cmd.getOptionValue("exp");
            }

            String geneannotation = null;
            if (cmd.hasOption("annotation")) {
                geneannotation = cmd.getOptionValue("annotation");
            }
            String output = null;
            if (cmd.hasOption("out")) {
                output = cmd.getOptionValue("out");
            }

            String input = null;
            if (cmd.hasOption("input")) {
                input = cmd.getOptionValue("input");
            }

            System.out.println(mode);
            switch (mode) {
                case "concatconditional":
                    Integer nriters = null;
                    if (cmd.hasOption("nriters")) {
                        nriters = (Integer) Integer.parseInt(cmd.getOptionValue("nriters"));
                    }
                    if (input == null || vcf == null || linkfile == null || nriters == null || output == null) {
                        System.out.println("Usage: --input qtlfilestr --vcf vcffile.vcf.gz --out output --gte samplelinkfile.txt --nriters nriters ");
                        System.out.println("Input: " + input);
                        System.out.println("VCF: " + vcf);
                        System.out.println("Output: " + output);
                        System.out.println("SampleLink: " + linkfile);
                        System.out.println("Nr Iters: " + nriters);
                    } else {
                        ConcatenateConditionalResults c = new ConcatenateConditionalResults();
                        c.run(input, vcf, linkfile, nriters, output);
                    }
                    break;
                case "determineld":
                    System.out.println("Determine LD between SNPs in two QTL files");
                    String input2 = null;
                    if (cmd.hasOption("input2")) {
                        input2 = cmd.getOptionValue("input2");
                    }
                    if (input == null || input2 == null || output == null || vcf == null) {
                        System.out.println("Usage: --input qtlfile1 --input2 qtlfile2 --vcf vcffile.vcf.gz --out output [--gte samplelinkfile.txt]");
                        System.out.println("Input: " + input);
                        System.out.println("Input2 " + input2);
                        System.out.println("Output: " + output);
                        System.out.println("VCF: " + vcf);
                        System.out.println("GTE: " + linkfile);
                    } else {
                        LDCalculator calc = new LDCalculator();
                        calc.run(input, input2, vcf, linkfile, output);
                    }
                    break;
                case "determineldgwas":
                    System.out.println("Determine LD between QTL SNPs and GWAS SNPs");
                    String eqtlset = null;
                    String eqtlfile = null;
                    String gwasset = null;
                    String gwasassocfile = null;
                    String gwaslistfile = null;
                    int prunedistance = 1000000;
                    double prunethreshold = 0.2;
                    double ldthreshold = 0.8;
                    boolean skipchr6 = false;
                    boolean matchbyrsid = false;

                    if (cmd.hasOption("eqtlset")) {
                        eqtlset = cmd.getOptionValue("eqtlset");
                    }
                    if (cmd.hasOption("eqtlfile")) {
                        eqtlfile = cmd.getOptionValue("eqtlfile");
                    }
                    if (cmd.hasOption("gwasset")) {
                        gwasset = cmd.getOptionValue("gwasset");
                    }
                    if (cmd.hasOption("gwasassocfile")) {
                        gwasassocfile = cmd.getOptionValue("gwasassocfile");
                    }
                    if (cmd.hasOption("gwaslistfile")) {
                        gwaslistfile = cmd.getOptionValue("gwaslistfile");
                    }
                    if (cmd.hasOption("prunedistance")) {
                        prunedistance = Integer.parseInt(cmd.getOptionValue("prunedistance"));
                    }
                    if (cmd.hasOption("prunethreshold")) {
                        prunethreshold = Double.parseDouble(cmd.getOptionValue("prunethreshold"));
                    }
                    if (cmd.hasOption("ldthreshold")) {
                        ldthreshold = Double.parseDouble(cmd.getOptionValue("ldthreshold"));
                    }
                    if (cmd.hasOption("skipchr6")) {
                        skipchr6 = true;
                    }
                    if (cmd.hasOption("matchbyrsid")) {
                        matchbyrsid = true;
                    }

                    if (eqtlset == null || eqtlfile == null
                            || gwasset == null || gwasassocfile == null || gwaslistfile == null
                            || vcf == null) {
                        System.out.println("Usage: --eqtlset eqtlset.txt.gz " +
                                "--eqtlfile eqtlfile.txt.gz " +
                                "--gwasset gwasset.txt.gz " +
                                "--gwasassocfile gwasassocfile.txt.gz " +
                                "--gwaslistfile gwaslist.txt.gz " +
                                "--vcf vcfCHR.vcf.gz" +
                                "[--gte samplelinkfile.txt] " +
                                "[--prunedistance 1E6] " +
                                "[--prunethreshold 0.2] " +
                                "[--ldthreshold 0.8] " +
                                "--out output.txt.gz " +
                                "[--skipchr6] " +
                                "[--matchbyrsid]");
                        System.out.println("eqtlset: " + eqtlset);
                        System.out.println("eqtlfile: " + eqtlfile);
                        System.out.println("gwasset: " + gwasset);
                        System.out.println("gwasassocfile: " + gwasassocfile);
                        System.out.println("gwaslistfile: " + gwaslistfile);
                        System.out.println("vcf: " + vcf);
                        System.out.println("gte: " + linkfile);
                        System.out.println("prunedistance: " + prunedistance);
                        System.out.println("prunethreshold: " + prunethreshold);
                        System.out.println("ldthreshold: " + ldthreshold);
                        System.out.println("out: " + output);
                        System.out.println("skipchr6: " + skipchr6);
                        System.out.println("matchbyrsid: " + matchbyrsid);

                    } else {
                        GWASLDCalculator calc = new GWASLDCalculator();

                        calc.runPerGWAS3(eqtlset, eqtlfile, gwasset, gwasassocfile, gwaslistfile, vcf, linkfile, prunedistance, prunethreshold, ldthreshold, output, skipchr6, matchbyrsid);

                    }
                    break;
                case "sortfile":
                    System.out.println("QTL file sorter sorts by position by default; use --sortz to sort by Z-score");
                    QTLFileSorter sorter = new QTLFileSorter();
                    if (input == null || output == null) {
                        System.out.println("Usage: --input infile.txt[.gz] --out outfile.txt[.gz] [--sortbyz|--sortbysnppos|--sortbygenepos]");
                        System.out.println("Input: " + input);
                        System.out.println("Output: " + output);
                        System.out.println("Sort by Z: " + cmd.hasOption("sortbyz"));
                    } else {
                        if (cmd.hasOption("sortbyz")) {
                            sorter.run(input, output, QTLFileSorter.SORTBY.Z);
                        } else if (cmd.hasOption("sortbysnppos")) {
                            sorter.run(input, output, QTLFileSorter.SORTBY.SNPPOS);
                        } else {
                            sorter.run(input, output, QTLFileSorter.SORTBY.GENEPOS);
                        }
                    }
                    break;
                case "regressqtl":
                    if (vcf == null || linkfile == null || geneannotation == null || genexpression == null || output == null || snpgenelimit == null) {
                        System.err.println("Usage: -m regressqtl --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ --snpgenelimit snpgenecombos.txt.gz");
                        System.err.println("Optional: --snpannotation snpannotation.txt.gz --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2 --genelimit");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Linkfile: " + linkfile);
                        System.out.println("Annotation: " + geneannotation);
                        System.out.println("Exp: " + genexpression);
                        System.out.println("Out: " + output);
                        System.out.println("SNP/Gene combos: " + snpgenelimit);
                        System.out.println("Genelimit: " + genelimit);
                    } else {

                        QTLRegression qtlr = new QTLRegression(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, minNrDatasets, minObservations, output);
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            qtlr.setReplaceMissingGenotypes(true);
                        }

                        if (cmd.hasOption("norank")) {
                            qtlr.setRankData(false);
                        }


                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            qtlr.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            qtlr.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            qtlr.setHwepthreshold(t);
                        }

                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            qtlr.setCisWindow(t);
                        }


                        if (cmd.hasOption("snpannotation")) {
                            qtlr.loadSNPAnnotation(cmd.getOptionValue("snpannotation"));
                        }
                        qtlr.run();
                    }
                    break;
                case "mbqtlplot":
                    if (vcf == null || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Required: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Chrom: " + chrom);
                        System.out.println("GTE:" + linkfile);
                        System.out.println("Gene annotation: " + geneannotation);
                        System.out.println("Gene expression: " + genexpression);
                        System.out.println("Output: " + output);
                    } else {
                        MbQTLPlot bpp = new MbQTLPlot(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, minNrDatasets, minObservations, output);
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            bpp.setReplaceMissingGenotypes(true);
                        }
                        if (cmd.hasOption("norank")) {
                            bpp.setRankData(false);
                        }
                        // bQTL.setMinObservations();
                        if (cmd.hasOption("outputall")) {
                            bpp.setOutputAll(true);
                        }


                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            bpp.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            bpp.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            bpp.setHwepthreshold(t);
                        }

                        bpp.plot();
                    }
                    break;
                case "mbqtl":
                    if (vcf == null || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Usage: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6 --nrdatasets 2");
                        System.err.println("Optional: --seed 123456789 --outputall --perm 1000 --snplog --outputallpermutations");
                        System.err.println("Optional: --expgroups expgroups.txt.gz");

                        System.out.println("You've set the following:");
                        System.out.println("VCF: " + vcf);
                        System.out.println("Chrom: " + chrom);
                        System.out.println("GTE:" + linkfile);
                        System.out.println("Gene annotation: " + geneannotation);
                        System.out.println("Gene expression: " + genexpression);
                        System.out.println("Output: " + output);
                    } else {
                        MbQTL2Parallel bQTL = new MbQTL2Parallel(vcf, chrom, linkfile, snplimit, genelimit, snpgenelimit, genexpression, geneannotation, minNrDatasets, minObservations, output);
                        if (cmd.hasOption("snpannotation")) {
                            bQTL.loadSNPAnnotation(cmd.getOptionValue("snpannotation"));
                        }

                        if (cmd.hasOption("testnonparseablechr")) {
                            bQTL.setTestNonParseableChr();
                        }

                        if (cmd.hasOption("expgroups")) {
                            bQTL.loadExpGroupAnnotation(cmd.getOptionValue("expgroups"));
                        }
                        if (cmd.hasOption("replacemissinggenotypes")) {
                            bQTL.setReplaceMissingGenotypes(true);
                        }
                        if (cmd.hasOption("norank")) {
                            bQTL.setRankData(false);
                        }
                        // bQTL.setMinObservations();
                        if (cmd.hasOption("outputall")) {
                            bQTL.setOutputAll(true);
                        }
                        if (cmd.hasOption("outputallpermutations")) {
                            bQTL.setOutputAllPermutations(true);
                        }
                        if (cmd.hasOption("snplog")) {
                            bQTL.setOutputSNPLog(true);
                        }

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            bQTL.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            bQTL.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            bQTL.setHwepthreshold(t);
                        }

                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            bQTL.setCisWindow(t);
                        }
                        if (cmd.hasOption("perm")) {
                            int t = Integer.parseInt(cmd.getOptionValue("perm"));
                            bQTL.setNrPermutations(t);
                        }
                        if (cmd.hasOption("seed")) {
                            long seed = Long.parseLong(cmd.getOptionValue("seed"));
                            bQTL.setRandomSeed(seed);
                        }


                        if (cmd.hasOption("empzmeta")) {
                            bQTL.setMetaanalysismethod(MetaAnalysisMethod.EMP);
                        }
                        if (cmd.hasOption("fisherzmeta")) {
                            bQTL.setMetaanalysismethod(MetaAnalysisMethod.FISHERZFIXED);
                        }
                        if (cmd.hasOption("fisherzmetarandom")) {
                            bQTL.setMetaanalysismethod(MetaAnalysisMethod.FISHERZRANDOM);
                        }

                        if (cmd.hasOption("trans")) {
                            bQTL.setAnalysisType(AnalysisType.TRANS);
                        } else if (cmd.hasOption("cistrans")) {
                            bQTL.setAnalysisType(AnalysisType.CISTRANS);
                        }

                        if (cmd.hasOption("mingenotypecount")) {
                            int mingenotypecount = Integer.parseInt(cmd.getOptionValue("mingenotypecount"));
                            bQTL.setMinGenotypeCount(mingenotypecount);
                        }

                        if (cmd.hasOption("splitmultiallelic")) {
                            bQTL.setSplitMultiAllelic();
                        }

                        if (cmd.hasOption("usehardgenotypecalls")) {
                            bQTL.setUseHardGenotypeCalls();
                        }

                        if (cmd.hasOption("onlytestsnps")) {
                            bQTL.setOnlyTestSNPs();
                        }

                        if (cmd.hasOption("correlationweights")) {
                            bQTL.setCorrelationWeights(cmd.getOptionValue("correlationweights"));
                        }

                        if (cmd.hasOption("permutationstrategy")) {
                            String setting = cmd.getOptionValue("permutationstrategy");
                            switch (setting) {
                                case "random":
                                    bQTL.setPermutationStrategy(PermutationStrategy.RANDOM);
                                    break;
                                case "semirandom":
                                    bQTL.setPermutationStrategy(PermutationStrategy.SEMIRANDOM);
                                    break;
                                case "matched":
                                    bQTL.setPermutationStrategy(PermutationStrategy.MATCHED);
                                    break;
                                default:
                                    System.err.println("Error: pick one permutation strategy: random, semimatched or matched");
                                    System.exit(-1);
                            }
                        }

                        bQTL.run();
                    }
                    break;
                case "mbqtlsingleds":
                    if (vcf == null || chrom == -1 || linkfile == null || geneannotation == null || genexpression == null || output == null) {
                        System.err.println("Usage: --vcf tabix.vcf.gz, --chr [1-22], --gte linkfile.txt, --annotation annotation.txt.gz, --exp expfile.txt.gz and --out /outdir/ ");
                        System.err.println("Optional: --replacemissinggenotypes, --norank, --minobservations 10 --maf 0.01 --cr 0.95 --hwep 0.001 --ciswindow 1E6");
                        System.err.println("Optional: --seed 123456789 --outputall --perm 1000 ");
                    } else {

                        MbQTLSingleDataset ds = new MbQTLSingleDataset(vcf, chrom, linkfile, genelimit, genexpression, geneannotation, minNrDatasets, minObservations, output);

                        if (cmd.hasOption("maf")) {
                            double t = Double.parseDouble(cmd.getOptionValue("maf"));
                            ds.setMafthreshold(t);
                        }
                        if (cmd.hasOption("cr")) {
                            double t = Double.parseDouble(cmd.getOptionValue("cr"));
                            ds.setCallratethreshold(t);
                        }
                        if (cmd.hasOption("hwep")) {
                            double t = Double.parseDouble(cmd.getOptionValue("hwep"));
                            ds.setHwepthreshold(t);
                        }
                        if (cmd.hasOption("ciswindow")) {
                            int t = Integer.parseInt(cmd.getOptionValue("ciswindow"));
                            ds.setCisWindow(t);
                        }
                        if (cmd.hasOption("perm")) {
                            int t = Integer.parseInt(cmd.getOptionValue("perm"));
                            ds.setNrPermutations(t);
                        }
                        ds.fastqtlclone2();
                    }
                    break;
                case "metaqtl":
                    System.out.println("Not implemented yet.");
                    break;
                default:

            }


        } catch (ParseException e) {
            System.out.println("Command line parse exception: " + e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.setWidth(160);
            formatter.printHelp("mbqtl.jar", options);
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
    }

}
