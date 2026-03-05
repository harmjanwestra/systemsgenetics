package mbqtl.gfx;

import java.util.Objects;

public class ForestPlotEQTL implements Comparable<ForestPlotEQTL> {
    public String EA;
    public String NONEA;
    public String snp;
    public boolean isMeta;
    String gene;
    public String dataset;
    public int samplesize;

    public double effectsize;
    public double error;

    public EFFECTTYPE effecttype = EFFECTTYPE.ZSCORE;

    public double pvalue = Double.NaN;

    public enum EFFECTTYPE {
        ZSCORE,
        CORRELATION,
        BETA
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ForestPlotEQTL eqtl = (ForestPlotEQTL) o;
        return Double.compare(eqtl.effectsize, effectsize) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(effectsize);
    }


    @Override
    public int compareTo(ForestPlotEQTL eqtl) {
        if (this.equals(eqtl)) {
            return 0;
        }

        if (this.effectsize > eqtl.effectsize) {
            return 1;
        } else {
            return -1;
        }


    }
}
