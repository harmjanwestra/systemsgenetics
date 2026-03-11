package mbqtl.stat;

import JSci.maths.statistics.NormalDistribution;
import cern.jet.stat.tdouble.Probability;
import umcg.genetica.math.stats.ZScores;

public class FisherWeightedMetaAnalysis {

    NormalDistribution d = new NormalDistribution();

    public double[] metaAnalyzeCorrelationsFixedEffect(double[] correlations, int[] sampleSizes) {
        return metaAnalyzeCorrelationsFixedEffect(correlations, sampleSizes, 0);
    }

    public double correlationToFisherZ(double r) {
        return 0.5 * Math.log((1 + r) / (1 - r)); // pp 42-43
    }

    /* From Borenstein et al, an introduction to Meta Analysis
    pages: 42, 66, 93
    */
    public double[] metaAnalyzeCorrelationsFixedEffect(double[] correlations, int[] sampleSizes, double Tsq) {

        double M = 0;
        double sWiYi = 0;
        double sWi = 0;
        double[] zvals = new double[correlations.length];
        for (int di = 0; di < correlations.length; di++) {
            double ri = correlations[di];
            if (!Double.isNaN(ri)) {

                double zi = 0.5 * Math.log((1 + ri) / (1 - ri)); // pp 42-43
                zvals[di] = zi;
                double ni = sampleSizes[di];
                double vi = 1d / (ni - 3);
                vi += Tsq; // Tsq is the
//                double se = Math.sqrt(vi);

                double wi = 1 / vi; // pp 66


                sWiYi += (wi * zi);
                sWi += wi;
            } else {
                zvals[di] = Double.NaN;
            }
        }

        M = sWiYi / sWi;
        double Vm = 1d / sWi;
        double SEm = Math.sqrt(Vm);

        double LLM = M - 1.96 * SEm;
        double ULM = M - 1.96 * SEm;

        double expr = Math.exp(2 * M);
        double metar = (expr - 1) / (expr + 1);
        double Z = M / SEm;
        double p = ZScores.zToP(Z);
        // 2 * (1 - d.cumulative(Math.abs(Z)));

        double[] het = determineHeterogeneity(zvals, sampleSizes, M);

        return new double[]{Z, p, metar, SEm, LLM, ULM, het[0], het[1], het[2], het[3]};

    }

    /* From Borenstein et al, an introduction to Meta Analysis
   pages: 145
   */
    public double[] determineHeterogeneity(double[] zvals, int[] sampleSizes, double M) {
        double Q = 0;
        int k = 0;
        for (int di = 0; di < zvals.length; di++) {
            double zi = zvals[di];
            if (!Double.isNaN(zi)) {
//                double zi = 0.5 * Math.log((1 + ri) / (1 - ri)); // pp 42-43
                double ni = sampleSizes[di];
                double vi = 1d / (ni - 3);
                double wi = 1d / vi;
                double zim = zi - M;
                Q += wi * zim * zim;
                k++;
            }
        }
        int df = k - 1;

        double ISq = (Q - df) / Q;
        ISq *= 100;
        double B;
        if (Q <= df + 1) {
            B = Math.sqrt(1d / (2 * (df - 1) * (1 - (1d / (3 * (df - 1) * (df - 1))))));
        } else {
            B = 0.5 * ((Math.log(Q) - Math.log(df)) / (Math.sqrt(2 * Q) - Math.sqrt(2 * df - 1)));
        }

        double L = Math.exp(0.5 * Math.log(Q / df) - 1.96 * B);
        double U = Math.exp(0.5 * Math.log(Q / df) + 1.96 * B);
        double LLi2 = (((L * L) - 1) / (L * L)) * 100; // Lower bound confidence interval
        double ULi2 = (((U * U) - 1) / (U * U)) * 100; // Upper bound confidence interval

        double p = Probability.chiSquareComplemented(df, Q);

        if (ISq < 0) {
            ISq = 0;
        }
        if (ISq > 100) {
            ISq = 100;
        }

        return new double[]{ISq, LLi2, ULi2, p};
    }

    /* From Borenstein et al, an introduction to Meta Analysis
    pages: 72
    */
    public double[] metaAnalyzeCorrelationsRandomEffect(double[] correlations, int[] sampleSizes) {

        double M = 0;
        double sWiYi = 0;
        double sWiYiSq = 0;
        int k = 0;
        double sWi = 0;
        double sWsq = 0;
        for (int di = 0; di < correlations.length; di++) {
            double ri = correlations[di];
            if (!Double.isNaN(ri)) {

                double zi = 0.5 * Math.log((1 + ri) / (1 - ri)); // pp 42-43
                double ni = sampleSizes[di];
                double vi = 1d / (ni - 3);
//                double se = Math.sqrt(vi);

                double wi = 1d / vi; // pp 66

                double zizi = zi * zi;
                sWiYi += (wi * zi);
                sWiYiSq += (wi * zizi);
                sWi += wi;
                sWsq += (wi * wi);
                k++;
            }
        }

        double Q = sWiYiSq - ((sWiYi * sWiYi) / sWi);
        double C = sWi - (sWsq / sWi);
        double df = k - 1;
        double Tsq = Q - (df / C);


        return metaAnalyzeCorrelationsFixedEffect(correlations, sampleSizes, Tsq);
    }


    // chapter 38 - Hunter and Schmidt

}
