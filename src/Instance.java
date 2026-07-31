import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Instance {
    public double gap = 0.00000000001;
    public double M = 1E10;
    public int v = 1;
    public int K = 6;
    public boolean isMovable = true;
    public int CGversion =2;
    public double emGAP = 0;
    public double C1rate = 0.2;
    public List<String[]> allNodeList;
    public double Tmax;
    public double W;
    public int numN;
    public int numC;
    public int numC1;
    public int numC2;
    public int[] C;
    public int[] C1;
    public int[] C2;
    public int[] Xcoord;
    public int[] Ycoord;
    public double[] q;
    public double[] P;
    public double totalP;
    public double ubP;
    public double delta = 0;
    public double[] d;
    public double[][] md;
    public int[][] feaMove;
    public double[][] rArc;
    public double[][] mArc;
    public double [][] T;
    public double [][] Trb;
    public double [][] Tmb;
    public double [] P_r;
    public double [] P_m;
    public double ps;
    public double[] mdMax;
    public double[] minMoveT;
    double total_Label = 0;
    int total_cgIter = 0;
    int total_NgS = 0;
    public long startT =  System.currentTimeMillis();
    public boolean overTime = false;
    public double[][] disMatrix;
    public int countPath = 0;
    
    public void readData(String filePath) {
        File file = null;
        FileInputStream fileInputStream = null;
        StringBuilder sb = null;
        try {
            file = new File(filePath);
            fileInputStream = new FileInputStream(file);

            int len;
            sb = new StringBuilder();
            while ((len = fileInputStream.read()) != -1) {
                sb.append((char) len);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (fileInputStream != null) {
                    fileInputStream.close();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        String[] data = sb.toString().split("\n");
        Pattern pattern = Pattern.compile("\\d+(\\.\\d+)?");
        Matcher matcher1 = pattern.matcher(data[data.length-1]);
        if (matcher1.find()){
            Tmax = Double.parseDouble(matcher1.group());
            System.out.println("Tmax="+Tmax);
        }
        allNodeList = new ArrayList<String[]>();
        for (int i = 1; i < data.length-1; i++) {
            String[] tempData = data[i].split("\\s+");
            allNodeList.add(tempData.clone());

        }
        allNodeList.add(data[1].split("\\s+"));
        numN = (int)allNodeList.size();
        numC = numN - 2;
        numC1 = (int) Math.ceil(numC * C1rate);
        numC2 = numC - numC1;
        C = java.util.stream.IntStream.rangeClosed(1, numN-2).toArray();
        C1 = java.util.stream.IntStream.rangeClosed(1, numC1).toArray();
        C2 = java.util.stream.IntStream.rangeClosed(numC1+1, numN-2).toArray();
        System.out.println("C="+ Arrays.toString(C));
        System.out.println("C1="+Arrays.toString(C1));
        System.out.println("C2="+Arrays.toString(C2));
        Xcoord = new int[numN];
        Ycoord = new int[numN];
        P = new double[numN];
        P_r = new double[numN];
        P_m = new double[numN];
        q = new double[numN];
        d = new double[numN];
        for (int i = 0; i < numN; i++) {
            Xcoord[i] = Integer.parseInt(allNodeList.get(i)[1]);
            Ycoord[i] = Integer.parseInt(allNodeList.get(i)[2]);
            q[i] = Double.parseDouble(allNodeList.get(i)[4]);
            P[i] = Double.parseDouble(allNodeList.get(i)[3]);
            P_r[i] = Double.parseDouble(allNodeList.get(i)[3]);
            P_m[i] = Double.parseDouble(allNodeList.get(i)[3]);
            d[i] = Double.parseDouble(allNodeList.get(i)[5]);
        }
        for (double qv : q) {
            W += qv;
        }
        W = Math.round(1.5*W/K);
        T = new double[numN][numN];
        for (int i = 0; i < numN; i++) {
            for (int j = i; j < numN; j++) {
                if (i != j) {
                    T[i][j] = Math.sqrt(Math.pow(Xcoord[i] - Xcoord[j],2) + Math.pow(Ycoord[i] - Ycoord[j],2)) / v;
                    T[i][j] = Math.round(T[i][j]*100) / 100.0;
                    T[j][i] = T[i][j];
                } else{
                    T[i][j] = M;
                }
            }
        }

        Trb = new double[numN][numN];
        Tmb = new double[numN][numN];
        for (int i = 0; i < numN; i++) {
            for (int j = 0; j < numN; j++) {
                Trb[i][j] = T[i][j];
                Tmb[i][j] = T[i][j];
            }
        }

        minMoveT = new double[numN];
        for (int i = 1; i < numN-1; i++) {
            double minValue = T[i][1];
            for (int j = 2; j < numN-1; j++) {
                if (T[i][j] < minValue) {
                    minValue = T[i][j];
                }
            }
            minMoveT[i] = minValue;
        }

        md = new double[numN][numN];
        if (isMovable) {
            for (int i = 1; i < numN-1; i++) {
                for (int j = 1; j < numN-1; j++) {

                    if (i != j) {
                        md[i][j] = Math.round((1 - T[i][j]/(T[0][i]+T[0][j]))*100)/100.0;

                        md[j][i] = md[i][j];
                    }
                }
            }
        } else {
            for (int i = 1; i < numN-1; i++) {
                for (int j = i; j < numN-1; j++) {
                    if (i != j) {
                        md[i][j] = Math.min(0,0);
                        md[j][i] = md[i][j];
                    }
                }
            }
        }

        mdMax = new double[numN];
        for (int i = 1; i < numN-1; i++) {
            double maxValue = md[i][1];
            for (int j = 2; j < numN-1; j++) {
                if (md[i][j] > maxValue) {
                    maxValue = md[i][j];
                }
            }
            mdMax[i] = maxValue;
        }

        feaMove = new int[numN][numN];
        for (int i = 1; i < numN-1; i++) {
            for (int j = 1; j < numN-1; j++) {
                if (i != j) {
                    if (P[j]*md[j][i] - T[j][i]*d[j] > 1e-7) {
                        feaMove[i][j] = 1;
                    } else if (P[j]*md[j][i] - T[j][i]*d[j] < -1e-7) {
                        feaMove[i][j] = 0;
                    }
                }
            }
        }

        rArc = new double[numN][numN];
        mArc = new double[numN][numN];
        for (int i = 0; i < P.length; i++) {
            totalP += P[i];
        }

        ubP = 0;
        for (int c : C) {
            double maxRP = P[c] - d[c]*T[0][c];
            double maxMP = 0;
            for (int m : C) {
                double tempMP = P[c]*md[c][m] - d[c]*Math.max(T[0][m], T[c][m]);
                if (tempMP > maxMP) {
                    maxMP = tempMP;
                }
            }
            ubP += Math.max(maxRP, maxMP);
        }

    }
    
    public void setPS(double ps) {
        this.ps = ps;
    }
}
