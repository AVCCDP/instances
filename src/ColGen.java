import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class ColGen {
    
    static class IloNumVarArray {
        int _num = 0;
        IloNumVar[] _array = new IloNumVar[32];
        void add(IloNumVar Var) {
            if (_num >= _array.length) {
                IloNumVar[] newArray = new IloNumVar[2 * _array.length];
                System.arraycopy(_array,0,newArray,0,_num);
                _array = newArray;
            }
            _array[_num++] = Var;
        }
        IloNumVar getVar(int i) {return _array[i];}
        int getSize() {
            return _num;
        }
    }

    
    public double runColGen(Instance it, ArrayList<scheme> schemes, int[] IntV) throws IOException {
        double obj = 0;
        double[] p1i;
        double[] p2i;
        boolean trigger;
        try {
            IloCplex MP = new IloCplex();
            IloObjective MPobj = MP.addMaximize();
            IloRange[] MPcst1 = new IloRange[it.numC1];
            for (int i = 0; i < it.numC1; i++) {
                MPcst1[i] = MP.addRange(1.0,1.0);
            }
            IloRange[] MPcst2 = new IloRange[it.numC2]
            for (int i : it.C2) {
                if (IntV[i] == 2) {
                    MPcst2[i-it.numC1-1] = MP.addRange(-Double.MAX_VALUE,1.0);
                }
                else if (IntV[i] == 1) {
                    MPcst2[i-it.numC1-1] = MP.addRange(1,1);
                }
                else if (IntV[i] == 0) {
                    MPcst2[i-it.numC1-1] = MP.addRange(0,0);
                }
            }

            IloRange MPcst3 = MP.addRange(0,it.K);
            IloNumVarArray MP_s = new IloNumVarArray();

            for (scheme s : schemes) {
                IloColumn cols = MP.column(MPobj, s.getProfit(it));
                int cstInd;
                for (Integer i : s.customers) {
                    if (i >= 1 && i < it.C2[0]) {
                        cstInd = i-1;
                        cols = cols.and(MP.column(MPcst1[cstInd],1));
                    }
                    else {
                        cstInd = i-it.numC1-1;
                        cols = cols.and(MP.column(MPcst2[cstInd],1));
                    }
                }

                cols = cols.and(MP.column(MPcst3,1.0));
                MP_s.add(MP.numVar(cols,0.0, Double.MAX_VALUE));
            }


            MP.setParam(IloCplex.IntParam.RootAlgorithm, IloCplex.Algorithm.Primal);
            MP.setOut(null);
            double[] allObj = new double[500];
            int iter = -1;

            DecimalFormat df = new DecimalFormat("#0000.00");
            trigger = true;
            boolean isEnum = false;
            while (trigger) {
                trigger = false;
                if (!MP.solve()) {
                    System.out.println("CG: MP infeasible!");
                    return 1E10;
                }
                allObj[++iter] = MP.getObjValue();
                System.out.print("\nCG Iter" + iter + " Current cost: "
                        + df.format(allObj[iter]) + " genSol:" + schemes.size());

                p1i = MP.getDuals(MPcst1);
                p2i = MP.getDuals(MPcst2);
                it.setPS(MP.getDual(MPcst3));
                for (int i : it.C1) {
                    it.P_r[i] = p1i[i-1] - it.P[i];
                    it.P_m[i] = p1i[i-1];
                }
                for (int i : it.C2) {
                    it.P_r[i] = p2i[i-it.numC1-1] - it.P[i];
                    it.P_m[i] = p2i[i-it.numC1-1];
                }

                ArrayList<scheme> newSchemes = new ArrayList<scheme>();
                System.out.println();
                System.out.println("查看gap="+ (it.ubP - MP.getObjValue())/it.ubP );
                ReLabeling sp = new ReLabeling();
                if (it.CGversion == 1) {
                    sp.GenMove(it);
                    sp.GenRouteV0(it,it.numC);
                    sp.GenScheme(newSchemes);
                    it.total_NgS += newSchemes.size();
                }
                else if (it.CGversion == 2) {
                    sp.GenMove(it);
                    sp.GenRouteV3(it,it.numC);
                    sp.GenScheme(newSchemes);
                    it.total_NgS += newSchemes.size();
                    if (newSchemes.isEmpty()) {
                        sp.GenMove(it);
                        sp.GenRouteV1(it,it.numC);
                        sp.GenScheme(newSchemes);
                        it.total_NgS += newSchemes.size();
                    }
                }
                sp = null;
                if (!newSchemes.isEmpty()) {
                    for (scheme s : newSchemes) {
                        IloColumn cols = MP.column(MPobj, s.getProfit(it));
                        int cstInd;
                        for (Integer i : s.customers) {
                            if (i >= 1 && i < it.C2[0]) {
                                cstInd = i-1;
                                cols = cols.and(MP.column(MPcst1[cstInd],1));
                            }
                            else {
                                cstInd = i-it.numC1-1;
                                cols = cols.and(MP.column(MPcst2[cstInd],1));
                            }
                        }
                        cols = cols.and(MP.column(MPcst3,1.0));
                        MP_s.add(MP.numVar(cols,0.0, Double.MAX_VALUE));
                        schemes.add(s);
                        trigger = true;
                    }
                    newSchemes = null;
                }
            }
            it.total_cgIter += (iter+1);

            for (int i = 0; i < MP_s.getSize(); i++) {
                schemes.get(i).setV(MP.getValue(MP_s.getVar(i)));
                if (Math.abs(schemes.get(i).v - 1) < 1E-7) {
                    schemes.get(i).setV(1);
                }
            }
            obj = MP.getObjValue();
            MP.end();
            return obj;

        } catch (IloException e) {
            System.err.println("Concert exception caught '" + e + "' caught");
        } catch (CloneNotSupportedException e) {
            throw new RuntimeException(e);
        }
        return 1E10;
    }
}
