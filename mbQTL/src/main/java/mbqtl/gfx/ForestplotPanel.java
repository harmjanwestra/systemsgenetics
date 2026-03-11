package mbqtl.gfx;

import umcg.genetica.graphics.DefaultGraphics;
import umcg.genetica.graphics.Range;
import umcg.genetica.graphics.panels.Panel;
import umcg.genetica.graphics.themes.DefaultTheme;

import java.awt.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ForestplotPanel extends Panel {

    private final int marginBetweenDots;
    private final int marginBetweenBoxes;
    private final int marginWithinBoxes;
    private final int dotsize;
    private final int minDotSize;
    private final int halfdotsize;
    private final int boxSize;
    private HashMap<String, ArrayList<ForestPlotEQTL>> data;
    private Color[] greyScale;
    private Color[] colors;
    private DefaultTheme deftheme;
    private ORIENTATION orientation = ORIENTATION.VERTICAL;
    private STATISTIC statistic = STATISTIC.ZSCORE;
    private Range rangeBeta;
    private boolean determineRangeOverAllEvents;
    private Graphics2D g2d;
    private boolean forceRange = false;
    private boolean shading = false;


    public void setStatistic(STATISTIC statistic){
        this.statistic = statistic;
    }

    public void setRange(Range range) {
        this.rangeBeta = range;
        forceRange = true;
    }

    public void setShading(boolean b) {
        shading = b;
    }

    private enum ORIENTATION {
        HORIZONTAL,
        VERTICAL
    }

    public enum STATISTIC {
        ZSCORE,
        CORRELATION,
        BETA
    }

    public void determineRangeOverAllEvents() {
        this.determineRangeOverAllEvents = true;
    }

    public ForestplotPanel(int nrRows, int nrCols) {
        super(nrRows, nrCols);

        marginBetweenDots = 10;
        marginBetweenBoxes = 20;
        marginWithinBoxes = 10;
        dotsize = 12;
        minDotSize = 4;
        halfdotsize = dotsize / 2;
        boxSize = 100;

        // define colors
        // metabrain color palette
        colors = new Color[]{
                new Color(0, 0, 0),
                new Color(230, 159, 0),
                new Color(86, 180, 233),
                new Color(204, 121, 167),
                new Color(240, 228, 66),
                new Color(0, 114, 178),
                new Color(213, 94, 0),
                new Color(0, 158, 115)
        };
        greyScale = new Color[]{
                new Color(0, 0, 0),
                new Color(128, 128, 128),
                new Color(169, 169, 169),
                new Color(211, 211, 211),
                new Color(232, 232, 232)
        };
    }

    public void setOrientation(ORIENTATION orientation) {
        this.orientation = orientation;
    }

    public void setData(ArrayList<ForestPlotEQTL> data) {
        HashMap<String, ArrayList<ForestPlotEQTL>> tmp = new HashMap<>();
        tmp.put(data.get(0).snp, data);
        this.data = tmp;
    }

    @Override
    public void draw(DefaultGraphics defaultGraphics) {
        g2d = defaultGraphics.getG2d();
        g2d.setColor(greyScale[0]);

        deftheme = (DefaultTheme) theme;
        deftheme.setColors(colors);

        ArrayList<String> keys = new ArrayList<>();
        keys.addAll(data.keySet());
        Collections.sort(keys);

        if (determineRangeOverAllEvents) {
            ArrayList<ForestPlotEQTL> tmp = new ArrayList<>();
            for (int i = 0; i < keys.size(); i++) {
                String outcome = keys.get(i);
                ArrayList<ForestPlotEQTL> events = data.get(outcome);
                tmp.addAll(events);
            }
            determineRanges(tmp);
        }

        //g2d.setFont(deftheme.getSmallFontBold());
        //            g2d.drawString("Gene", x0 - getStringWidth(g2d, "Gene") - marginBetweenDots, y1 - 3);

        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[1]);
            g2d.setFont(deftheme.getSmallFont());
            if(statistic==STATISTIC.BETA){
                int strw = getStringWidth(g2d, "Beta + 95% CI");
                 strw = strw/2;
                g2d.drawString("Beta + 95% CI", x0 + (boxSize / 2) - strw, y0 - 20);
            } else if(statistic==STATISTIC.CORRELATION){
                int strw = getStringWidth(g2d, "Correlation + s.e.");
                strw = strw/2;
                g2d.drawString("Correlation + s.e.", x0 + (boxSize / 2) - strw, y0 - 20);
            } else {
                int strw = getStringWidth(g2d, "Z-score");
                strw = strw/2;
                g2d.drawString("Z-score", x0 + (boxSize / 2) - strw, y0 - 20);
            }

        }

        int startX = x0;
        int startY = y0;
        int octr = 0;
        for (int i = 0; i < keys.size(); i++) {
            String outcome = keys.get(i);
            ArrayList<ForestPlotEQTL> events = data.get(outcome);
            if (!determineRangeOverAllEvents) {
                rangeBeta = null;
                determineRanges(events);
            }

            System.out.println("Range Beta: " + rangeBeta);

            plotData(events, startX, startY, outcome, defaultGraphics);
            if (orientation == ORIENTATION.HORIZONTAL) {
                int boxWidth = events.size() * (dotsize + marginBetweenDots);
                startX += boxWidth + marginBetweenBoxes;
            } else {
                int boxHeight = events.size() * (dotsize + marginBetweenDots) + (2 * marginWithinBoxes);
                startY += boxHeight + marginBetweenBoxes;
            }
            octr++;

        }
    }

    private void drawBox(int x1, int y1, int boxHeight, int boxWidth) {

        // draw 0 line
        g2d.setColor(greyScale[2]);
        Stroke defStroke = g2d.getStroke();
        g2d.setStroke(deftheme.strokeDashed);
        if (orientation == ORIENTATION.HORIZONTAL) {
            g2d.drawLine(x1, y1 + getPos(0, boxHeight, rangeBeta), x1 + boxWidth, y1 + getPos(0, boxHeight, rangeBeta));
        } else {
            g2d.drawLine(x1 + getPos(0, boxWidth, rangeBeta), y1, x1 + getPos(0, boxWidth, rangeBeta), y1 + boxHeight);
        }
        g2d.setStroke(defStroke);

        g2d.setColor(greyScale[0]);
        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[2]);
        }
        g2d.drawLine(x1, y1, x1, y1 + boxHeight); // |--
        g2d.drawLine(x1 + boxWidth, y1, x1 + boxWidth, y1 + boxHeight); // --|

        g2d.setColor(greyScale[4]);
        if (orientation == ORIENTATION.VERTICAL) {
            g2d.setColor(greyScale[0]);
        }
        g2d.drawLine(x1, y1, x1 + boxWidth, y1); // --
        g2d.drawLine(x1, y1 + boxHeight, x1 + boxWidth, y1 + boxHeight); // __


    }

    private void plotData(ArrayList<ForestPlotEQTL> events, int x1, int y1, String outcome, DefaultGraphics defaultGraphics) {
        int boxWidth = 0;
        int boxHeight = 0;
        if (orientation == ORIENTATION.HORIZONTAL) {
            boxHeight = boxSize;
            boxWidth = events.size() * (dotsize + marginBetweenDots);
            drawBox(x1, y1, boxHeight, boxWidth); // draw beta box

            int y2 = y1 + boxHeight + marginBetweenBoxes;
            drawBox(x1, y2, boxHeight, boxWidth); // draw WR box


            // plot beta values
            int ectr = 0;
            for (ForestPlotEQTL event : events) {
                g2d.setColor(colors[1]);
                // plot Beta
                int posX = ectr * (dotsize + marginBetweenDots);
                g2d.fillOval(x1 + posX, y1 + getPos(event.effectsize, boxHeight, rangeBeta), dotsize, dotsize);


                // event snp
                // event gene
                // event allele
            }
        } else {
            boxWidth = boxSize;
            boxHeight = events.size() * (dotsize + marginBetweenDots) + (2 * marginWithinBoxes);

            if (shading) {
                // draw shading every 5 lines
                if (events.size() > 10) {
                    boolean shade = true;
                    for (int i = 0; i < events.size(); i++) {
                        if (i > 0 && i % 5 == 0) {
                            if (shade) {
                                // shade the next 5 genes
                                int stop = i + 5;
                                if (stop > events.size()) {
                                    stop = events.size();
                                }
                                int posY = y1 + i * (dotsize + marginBetweenDots) + (dotsize / 2) + marginWithinBoxes;
                                int posY2 = y1 + stop * (dotsize + marginBetweenDots) + (dotsize / 2) + marginWithinBoxes;
                                g2d.setColor(greyScale[4]);
                                g2d.fillRect(0, posY, width, posY2 - posY);
                                shade = false;
                            } else {
                                shade = true;
                            }
                        }
                    }
                }
            }

            // plot the scales
            g2d.setColor(greyScale[1]);
            DecimalFormat df = new DecimalFormat("#.##");
            g2d.drawString(df.format(rangeBeta.getMinX()), x1 - (getStringWidth(g2d, df.format(rangeBeta.getMinX())) / 2), y1 - 5);
            g2d.drawString(df.format(0), x1 + getPos(0, boxWidth, rangeBeta) - (getStringWidth(g2d, "0") / 2), y1 - 5);
            g2d.drawString(df.format(rangeBeta.getMaxX()), x1 + boxWidth - (getStringWidth(g2d, df.format(rangeBeta.getMaxX())) / 2), y1 - 5);


            drawBox(x1, y1, boxHeight, boxWidth); // draw beta box
            int x2 = x1;



            g2d.setFont(deftheme.getSmallFont());
            int smallFontHeight = g2d.getFontMetrics().getHeight();

            int maxSNPStrWidth = 0;
            int maxGeneStrWidth = 0;
            for (ForestPlotEQTL event : events) {
                int strWidth = super.getStringWidth(g2d, event.EA + " - " + event.snp);
                if (strWidth > maxSNPStrWidth) {
                    maxSNPStrWidth = strWidth;
                }
                if (event.gene != null) {
                    int strWidthGene = super.getStringWidth(g2d, event.gene);
                    if (strWidthGene > maxGeneStrWidth) {
                        maxGeneStrWidth = strWidthGene;
                    }
                }
            }

            // plot outcome name
            g2d.setColor(greyScale[0]);
            g2d.setFont(deftheme.getSmallFontBold());
            int x3 = x1 + boxWidth + marginBetweenDots;

            g2d.drawString(outcome, x3, y1 - 3 + marginWithinBoxes);

            g2d.setFont(deftheme.getSmallFont());
            // plot beta values
            int ectr = 0;

            int maxN = 0;
            for (ForestPlotEQTL event : events) {
                if (!event.isMeta) {
                    if (event.samplesize > maxN) {
                        maxN = event.samplesize;
                    }
                }
            }

            for (ForestPlotEQTL event : events) {

                if (event.isMeta) {
                    g2d.setColor(colors[2]);

                } else {
                    g2d.setColor(colors[1]);
                }

                int dotSizeForEvent = dotsize;
                if (!event.isMeta) {
                    dotSizeForEvent = (int) Math.ceil(dotsize * ((double) event.samplesize / maxN));
                    if (dotSizeForEvent < minDotSize) {
                        dotSizeForEvent = minDotSize;
                    }

                }
                // plot Beta exposure
                int posY = ectr * (dotsize + marginBetweenDots) + dotsize + marginWithinBoxes;
                if (event.isMeta) {
                    int xp0 = x1 + getPos(event.effectsize, boxWidth, rangeBeta) - (dotsize / 2);
                    int yp0 = y1 + posY + (dotsize / 2);
                    Polygon p = new Polygon();
                    p.addPoint(xp0, yp0);
                    p.addPoint(xp0 + (dotsize / 2), yp0 - (dotsize / 2));
                    p.addPoint(xp0 + dotsize, yp0);
                    p.addPoint(xp0 + (dotsize / 2), yp0 + (dotsize / 2));
                    g2d.fillPolygon(p);
                } else {
                    g2d.fillOval(x1 + getPos(event.effectsize, boxWidth, rangeBeta) - (dotSizeForEvent / 2), y1 + posY, dotSizeForEvent, dotSizeForEvent);
                }


                // plot beta exposure confidence
                double conf1 = event.effectsize;// - event.error; // confinterval = event.effectsize - 1.96 * event.error
                double conf2 = event.effectsize;// + event.error;
                if(statistic == STATISTIC.BETA){
                    conf1 = event.effectsize - 1.96 * event.error; // confinterval = event.effectsize - 1.96 * event.error
                    conf2 = event.effectsize + 1.96 * event.error;
                } else if (statistic == STATISTIC.CORRELATION){
                    conf1 = event.effectsize - event.error; // confinterval = event.effectsize - 1.96 * event.error
                    conf2 = event.effectsize + event.error;
                }

                g2d.drawLine(x1 + getPos(conf1, boxWidth, rangeBeta), y1 + posY + (dotSizeForEvent / 2), x1 + getPos(conf2, boxWidth, rangeBeta), y1 + posY + (dotSizeForEvent / 2));

                // event snp
                g2d.setColor(greyScale[1]);
                g2d.setFont(deftheme.getSmallFont());

                // dataset name
                posY = ectr * (dotsize + marginBetweenDots) + dotsize + marginBetweenDots + marginWithinBoxes; // + (smallFontHeight);

                DecimalFormat f = new DecimalFormat("#.###");
                DecimalFormat pf = new DecimalFormat("0.0E0");
                String dsStr = event.dataset;
                String pstr = "";
                String zstr = "";

                double p = event.pvalue;
                if (!Double.isNaN(p)) {
                    if (p > 0.01) {
                        pstr = f.format(p);
                    } else {
                        pstr = pf.format(p);
                    }
                    pstr="; p="+pstr;
                }
                if(!Double.isNaN(event.effectsize)){
                    if(statistic == STATISTIC.BETA){
                        zstr = "; b="+f.format(event.effectsize) +" se="+f.format(event.error);
                    } else if(statistic == STATISTIC.CORRELATION){
                        zstr = "; r="+f.format(event.effectsize) +" se="+f.format(event.error);
                    } else {
                        zstr = "; z="+f.format(event.effectsize);
                    }
                }

                g2d.drawString(event.dataset + "\t(n="+event.samplesize+zstr+pstr+")", x3, y1 + posY);


                // event gene
                if (event.gene != null) {
                    g2d.drawString(event.gene, x0 - getStringWidth(g2d, event.gene) - marginBetweenDots, y1 + posY);
                }
                ectr++;

            }
        }
    }

    private int getPos(double val, int size, Range range) {
        double percB = range.getRelativePositionX(val);
        int pixels = (int) Math.ceil(percB * size);
        return pixels;
    }

    private void determineRanges(ArrayList<ForestPlotEQTL> events) {
        if (!forceRange) {
            // determine min, max
            double minBeta = Double.MAX_VALUE;
            double maxBeta = -Double.MAX_VALUE;
            for (ForestPlotEQTL e : events) {

                double valmin = e.effectsize - (1.96 * e.error);
                double valmax = e.effectsize + (1.96 * e.error);

                if (valmin < minBeta) {
                    minBeta = valmin;
                }
                if (valmax > maxBeta) {
                    maxBeta = valmax;
                }
            }

            // make symmetrical
            minBeta = -Math.max(Math.abs(minBeta), Math.abs(maxBeta));
            maxBeta = -minBeta;

            // determine ranges of values and round
            rangeBeta = new Range(minBeta, minBeta, maxBeta, maxBeta);
            rangeBeta.roundX();
            rangeBeta.roundY();
        }

    }

    public void setData(HashMap<String, ArrayList<ForestPlotEQTL>> eqtls) {
        this.data = eqtls;
    }
}
