package mbqtl;

public class MbQTLDebug {

    public static void main(String[] args) {

        String[] args2 = new String[]{
                "--mode", "mbqtl",
                "--vcf", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/genotypes/chrCHR.vcf.gz",
                "--gte", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/gte.txt",
                "--exp", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/brain-exp.txt.gz",
                "--annotation", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/ProbeAnnotation-gencode.v45.primary_assembly.annotation.txt.gz",
                "--out", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/output/test",
                "--cis",
                "--perm", "0"
        };

//        "--snpgenelimit","/Users/harm-jan/SyncData/TMP/orfeas/transtest/snpGeneCombos.txt.gz",
//                "--snpannotation", "/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/genotypes/snpannotation.txt.gz",
        // "--snpannotation","/Users/harm-jan/SyncData/TMP/orfeas/transtest/data/genotypes/chrCHR.vcf.gz",


        Main.main(args2);


    }
}
