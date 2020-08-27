package org.pankratzlab.internal.gwas;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
// import java.io.*;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.gui.UITools;
import org.pankratzlab.common.mining.Transformations;
import org.pankratzlab.utils.Grafik;

public class MatchesVisualized {

  public static final int WIDTH_BUFFER = 25;
  public static final int HEIGHT_BUFFER = 25;
  public static final int SIZE = 8;
  private final String[] anchors;
  private final String[] barnacles;
  private double[][] data;
  private final int[][] pairs;
  private final double[] dists;
  private final int x = 0;
  private final int y = 1;

  public MatchesVisualized(String dir, String anchorList, String barnacleList, String factorfile,
                           int[] factorIndices, String pairings) {
    String[] line;
    Hashtable<String, String> hash;
    Vector<String> v;
    long time;
    double[][] trans;

    time = new Date().getTime();
    anchors = HashVec.loadFileToStringArray(dir + anchorList, false, new int[] {0}, true);
    barnacles = HashVec.loadFileToStringArray(dir + barnacleList, false, new int[] {0}, true);
    hash = HashVec.loadFileToHashString(dir + factorfile, 0, factorIndices, "\t", true);

    data = new double[anchors.length + barnacles.length][factorIndices.length];
    for (int i = 0; i < anchors.length; i++) {
      data[i] = ArrayUtils.toDoubleArray(hash.get(anchors[i]).split(PSF.Regex.GREEDY_WHITESPACE));
    }
    for (int i = 0; i < barnacles.length; i++) {
      data[anchors.length + i] = ArrayUtils.toDoubleArray(hash.get(barnacles[i])
                                                              .split(PSF.Regex.GREEDY_WHITESPACE));
    }
    trans = Matrix.transpose(data);
    for (int i = 0; i < factorIndices.length; i++) {
      trans[i] = Transformations.standardizeRange(trans[i]);
    }
    data = Matrix.transpose(trans);

    v = HashVec.loadFileToVec(dir + pairings, true, false, false);
    if (v.size() != anchors.length) {
      System.err.println("Error - number of pairings (" + v.size()
                         + ") doesn't match number of anchors loaded (" + anchors.length + ")");
      System.exit(1);
    }

    pairs = new int[anchors.length][2];
    dists = new double[anchors.length];
    for (int i = 0; i < pairs.length; i++) {
      line = v.elementAt(i).split(PSF.Regex.GREEDY_WHITESPACE);
      pairs[i][0] = ext.indexOfStr(line[0], anchors);
      pairs[i][1] = ext.indexOfStr(line[1], barnacles);
      try {
        dists[i] = Double.parseDouble(line[2]);
      } catch (NumberFormatException nfe) {
        dists[i] = Double.NaN;
      }

    }

    JFrame frame = new JFrame(ext.rootOf(pairings));
    frame.setMinimumSize(new Dimension(20, 20));
    UITools.setSize(frame, new Dimension(1000, 720));
    frame.setVisible(true);

    JPanel panel = new JPanel() {

      public static final long serialVersionUID = 7L;

      @Override
      public void paintComponent(Graphics g) {
        double mean = ArrayUtils.mean(dists, true);
        double stdev = ArrayUtils.stdev(dists);

        for (int i = 0; i < pairs.length; i++) {
          if (!Double.isNaN(dists[i])) {
            Grafik.drawThickLine(g,
                                 (int) (data[pairs[i][0]][x] * (getWidth() - 2 * WIDTH_BUFFER))
                                    + WIDTH_BUFFER,
                                 getHeight() - (int) (data[pairs[i][0]][y]
                                                      * (getHeight() - 2 * HEIGHT_BUFFER))
                                                    - HEIGHT_BUFFER,
                                 (int) (data[anchors.length + pairs[i][1]][x]
                                        * (getWidth() - 2 * WIDTH_BUFFER)) + WIDTH_BUFFER,
                                 getHeight() - (int) (data[anchors.length
                                                           + pairs[i][1]][y]
                                                      * (getHeight() - 2 * HEIGHT_BUFFER))
                                                                                           - HEIGHT_BUFFER,
                                 4, dists[i] < mean + 3 * stdev ? Color.BLUE : Color.ORANGE);
          }
        }

        g.setColor(Color.RED);
        for (int i = 0; i < anchors.length; i++) {
          g.fillOval((int) (data[i][x] * (getWidth() - 2 * WIDTH_BUFFER)) + WIDTH_BUFFER - SIZE / 2,
                     getHeight() - (int) (data[i][y] * (getHeight()
                                                        - 2 * HEIGHT_BUFFER)) - HEIGHT_BUFFER - SIZE
                                                                                                / 2,
                     SIZE, SIZE);
        }
        g.setColor(Color.BLACK);
        for (int i = anchors.length; i < data.length; i++) {
          g.fillOval((int) (data[i][x] * (getWidth() - 2 * WIDTH_BUFFER)) + WIDTH_BUFFER - SIZE / 2,
                     getHeight() - (int) (data[i][y] * (getHeight()
                                                        - 2 * HEIGHT_BUFFER)) - HEIGHT_BUFFER - SIZE
                                                                                                / 2,
                     SIZE, SIZE);
        }
      }
    };
    frame.getContentPane().add(panel);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.pack();

    System.out.println("Finished writing distances_" + ArrayUtils.toStr(factorIndices, ",") + " in "
                       + ext.getTimeElapsed(time));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\gwas\\MatchingForMito\\";
    String anchors = "anchor_cases.dat";
    String barnacles = "barnacle_controls.dat";
    String factors = "mds10.mds.xln";
    // String pairings = "dsts_norm_minMin.xln";
    // String pairings = "dsts_norm_maxMin.xln";
    // String pairings = "distances_1,2_norm_minMin.xln";
    // String pairings = "distances_1,2_norm_maxMin.xln";
    // String pairings = "distances_1-100_norm_maxMin.xln";
    // String pairings = "distances_C1_normx4,C2_normx4,Age_normx4,Malex1_norm_minMin.xln";
    String pairings = "distances_C1_normx8,C2_normx8,Age_normx4,Malex1_norm_minMin.xln";

    // int[] factorIndices = new int[] {1,2,3,4};
    // int[] factorIndices = new int[] {1,2,3,4,5,6,7,8,9,10};
    int[] factorIndices = new int[] {1, 2};

    String usage = "\\n" + "kaput.MatchesVisualized requires 0-1 arguments\n"
                   + "   (0) directory (i.e. dir=" + dir + " (default))\n"
                   + "   (1) anchors (i.e. anchors=" + anchors + " (default))\n"
                   + "   (2) barnacles (i.e. barnacles=" + barnacles + " (default))\n"
                   + "   (3) file with factors (i.e. factors=" + factors + " (default))\n"
                   + "   (4) indices of factors in clusterfile (i.e. indices="
                   + ArrayUtils.toStr(factorIndices, ",") + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("anchors=")) {
        anchors = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("barnacles=")) {
        barnacles = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("factors=")) {
        factors = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("indices=")) {
        factorIndices = ArrayUtils.toIntArray(arg.split("=")[1].split(","));
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new MatchesVisualized(dir, anchors, barnacles, factors, factorIndices, pairings);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}