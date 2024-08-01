package mbqtl.datastructures;

import umcg.genetica.enums.Chromosome;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.text.Strings;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class SNPAnnotation {
    HashMap<String, Integer> variantToId = new HashMap<>();
    ArrayList<String> variants = new ArrayList<>();
    ArrayList<Chromosome> chromosomes = new ArrayList<>();
    ArrayList<Integer> positions = new ArrayList<>();

    Set<String> snpLimit = null;

    public SNPAnnotation(String snpannotation) throws IOException {
        load(snpannotation);
    }

    public SNPAnnotation(String snpannotation, Set<String> snpLimit) throws IOException {
        this.snpLimit = snpLimit;
        load(snpannotation);
    }

    private void load(String snpannotation) throws IOException {
        System.out.println("Loading variant annotation from: " + snpannotation);

        if (snpannotation.endsWith(".vcf.gz") || snpannotation.endsWith(".vcf")) {
            if (snpannotation.contains("CHR")) {
                System.out.println("SNP annotation is VCF with CHR template string.");
                for (int chr = 1; chr < 23; chr++) {
                    String tmpVCF = snpannotation.replace("CHR", "" + chr);
                    if (Gpio.exists(tmpVCF)) {
                        loadFromVCF(tmpVCF);
                    } else {
                        System.out.println("Notice: VCF file for chr" + chr + " not found.");
                    }
                }
            } else {
                if (Gpio.exists(snpannotation)) {
                    loadFromVCF(snpannotation);
                }
            }
            for (int v = 0; v < variants.size(); v++) {
                String name = variants.get(v);
                if (!variantToId.containsKey(name)) {
                    variantToId.put(variants.get(v), v);
                } else {
                    System.out.println("Duplicate variant ID in VCF: " + name);
                }
            }
        } else {
            if (Gpio.exists(snpannotation)) {
                loadFromTxt(snpannotation);
            } else {
                System.out.println("Provided SNP annotation file does not exist: " + snpannotation);
            }
        }

        if (snpLimit != null) {
            System.out.println("Loaded annotation for " + variants.size() + " variants out of " + snpLimit.size() + " in SNP limit set.");
        } else {
            System.out.println("Loaded annotation for " + variants.size() + " variants.");
        }

    }

    private void loadFromVCF(String vcf) throws IOException {
        System.out.println("Parsing: " + vcf);
        TextFile tf = new TextFile(vcf, TextFile.R);
        String ln = tf.readLine();
        while (ln != null) {
            String[] elems = Strings.subsplit(ln, Strings.tab, 0, 4);
            String variantid = elems[2];
            if (snpLimit == null || snpLimit.contains(variantid)) {
                int pos = Integer.parseInt(elems[1]);
                Chromosome chr = Chromosome.parseChr(elems[0]);
                variants.add(variantid);
                positions.add(pos);
                chromosomes.add(chr);
            }
            ln = tf.readLine();
        }
        tf.close();
    }

    private void loadFromTxt(String snpannotation) throws IOException {
        TextFile tf = new TextFile(snpannotation, TextFile.R);
        String[] elems = tf.readLineElems(TextFile.tab);
        int ctr = 0;
        while (elems != null) {
            if (elems.length > 2) {
                try {
                    String variant = elems[0];
                    if (snpLimit == null || snpLimit.contains(variant)) {
                        Chromosome chr = Chromosome.parseChr(elems[1]);
                        int pos = Integer.parseInt(elems[2]);

                        if (variantToId.containsKey(variant)) {
                            int tmpid = variantToId.get(variant);
                            Chromosome tmpchr = chromosomes.get(tmpid);
                            Integer tmppos = positions.get(tmpid);
                            if (!tmpchr.equals(chr) || !tmppos.equals(pos)) {
                                System.out.println("Variant " + variant + " has multiple annotations. Is this correct? Only storing first instance.");
                            }
                        } else {
                            variantToId.put(variant, ctr);
                            variants.add(variant);
                            chromosomes.add(chr);
                            positions.add(pos);
                            ctr++;
                        }
                    }

                } catch (NumberFormatException e) {

                }
            }
            elems = tf.readLineElems(TextFile.tab);
        }
        tf.close();

    }

    public Integer getId(String variantName) {
        return variantToId.get(variantName);
    }

    public int getPos(Integer id) {
        return positions.get(id);
    }

    public Chromosome getChr(Integer id) {
        return chromosomes.get(id);
    }
}
