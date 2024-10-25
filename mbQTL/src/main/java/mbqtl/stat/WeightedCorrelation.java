package mbqtl.stat;


public class WeightedCorrelation {

    private static double partialCov(double xi, double mwx, double yi, double mwy, double wi) {
        double xminxw = xi - mwx;
        double yminyw = yi - mwy;
        return xminxw * yminyw * wi;
    }

    public static double correlate(double[] x, double[] y, double[] w) {
        double sumw = 0;
        double mwx = 0;
        double mwy = 0;

        // weighted mean
        for (int i = 0; i < x.length; i++) {
            mwx += x[i] * w[i];
            mwy += y[i] * w[i];
            sumw += w[i];
        }

        mwx /= sumw;
        mwy /= sumw;


        // covariances
        double covxyw = 0;
        double covyyw = 0;
        double covxxw = 0;
        for (int i = 0; i < x.length; i++) {

            double xi = x[i];
            double yi = y[i];
            double wi = w[i];

            covxyw += partialCov(xi, mwx, yi, mwy, wi);
            covxxw += partialCov(xi, mwx, xi, mwx, wi);
            covyyw += partialCov(yi, mwy, yi, mwy, wi);

        }

        covxyw /= sumw;
        covxxw /= sumw;
        covyyw /= sumw;

        return covxyw / Math.sqrt(covxxw * covyyw);
    }

}
