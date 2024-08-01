package mbqtl.vcf;

import mbqtl.vcf.filter.genotypefilters.VCFGenotypeFilter;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Created by hwestra on 5/11/15.
 */
public class VCFGenotypeData implements Iterator<VCFVariant> {

    private Set<String> variantLimitSet;
    private ArrayList<VCFGenotypeFilter> genotypeFilters;
    ArrayList<String> samples = new ArrayList<String>();
    TextFile tf = null;
    String next = null;
    boolean[] samplefilter = null;
    HashSet<String> includeTheseSamples = null;
    private boolean hasnext = true;

    public VCFGenotypeData(TextFile tf, HashSet<String> includeTheseSamples, ArrayList<VCFGenotypeFilter> genotypeFilters) throws IOException {
        this.tf = tf;
        this.genotypeFilters = genotypeFilters;
        this.includeTheseSamples = includeTheseSamples;
        tf.close();
        tf.open();

        init();
    }

    public VCFGenotypeData(TextFile tf, HashSet<String> includeTheseSamples) throws IOException {
        this.tf = tf;
        this.includeTheseSamples = includeTheseSamples;
        tf.close();
        tf.open();

        init();
    }

    public VCFGenotypeData(TextFile tf, boolean[] sampleFilter) throws IOException {
        this.tf = tf;
        this.samplefilter = sampleFilter;
        tf.close();
        tf.open();

        init();
    }

    public VCFGenotypeData(TextFile tf) throws IOException {
        this.tf = tf;
        tf.close();
        tf.open();
        init();
    }

    public VCFGenotypeData(String vcffile) throws IOException {
        tf = new TextFile(vcffile, TextFile.R);
        init();

    }

    public VCFGenotypeData(String vcfFile, boolean[] genotypeSamplesToInclude, Set<String> variantLimitSet) throws IOException {
        this.samplefilter = genotypeSamplesToInclude;
        this.variantLimitSet = variantLimitSet;
        tf = new TextFile(vcfFile, TextFile.R);
        init();
    }

    public void open() throws IOException {
        tf.open();
        init();
    }

    private void init() throws IOException {
        String ln = tf.readLine();
        while (ln.startsWith("#")) {
            if (ln.startsWith("#CHROM")) {
                parseheader(ln);
                break;
            }
            ln = tf.readLine();
        }
    }


    public ArrayList<String> getSamples() {
        return samples;
    }

    private void parseheader(String ln) {
        String[] elems = ln.split("\t");
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

        for (int e = 9; e < elems.length; e++) {
            String sample = elems[e];
            samples.add(sample);


        }
    }

    @Override
    public boolean hasNext() {
        return hasnext;
    }

    @Override
    public VCFVariant next() {
        try {
            String next = tf.readLine();
            if (next == null) {
                hasnext = false;
                return null;
            } else {
                if (variantLimitSet != null) {
                    if (variantLimitSet.contains(new VCFVariant(next, VCFVariant.PARSE.HEADER, samplefilter).getId())) {
                        return new VCFVariant(next, VCFVariant.PARSE.ALL, samplefilter);
                    } else {
                        return null;
                    }
                } else {
                    return new VCFVariant(next, VCFVariant.PARSE.ALL, samplefilter);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        hasnext = false;
        return null;

//
//        if (next != null) {
//            if (variantLimitSet != null) {
//
//            }
//            VCFVariant var = new VCFVariant(next, VCFVariant.PARSE.ALL, samplefilter);
//            try {
//                next = tf.readLine();
//            } catch (IOException ex) {
//                throw new RuntimeException(ex);
//            }
//        }
//        return null;

    }

    @Override
    public void remove() {

    }

    public void close() throws IOException {
        this.tf.close();
    }

//	public VCFVariant nextLoadGenotypesOnly() {
//		VCFVariant current = next;
//		try {
//
//			String ln = tf.readLine();
//			if (ln != null) {
//				next = new VCFVariant(ln, VCFVariant.PARSE.GENOTYPES);
//			} else {
//				next = null;
//			}
//
//		} catch (IOException ex) {
//			throw new RuntimeException(ex);
//		}
//		return current;
//	}

    public VCFVariant nextLoadHeader() {
        VCFVariant current = new VCFVariant(next, VCFVariant.PARSE.HEADER);
        try {

            String ln = tf.readLine();
            if (ln != null) {
                next = ln; // new VCFVariant(ln, VCFVariant.PARSE.HEADER);
            } else {
                next = null;
            }

        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
//		return current;
        return null;
    }


}
