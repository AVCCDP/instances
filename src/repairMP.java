import ilog.concert.*;
import ilog.cplex.IloCplex;

import java.io.IOException;
import java.util.ArrayList;

public class repairMP {
    double IP_obj;
    double PS;
    
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

    public double runIP(Instance it, ArrayList<scheme> schemes, int[] IntV) throws IOException {
        double obj;
        try {
            IloCplex IP = new IloCplex();
            IloObjective IPobj = IP.addMaximize();
            IloRange[] IPcst1 = new IloRange[it.numC1];
            for (int i = 0; i < it.numC1; i++) {
                IPcst1[i] = IP.addRange(1.0,1.0);
            }
            IloRange[] IPcst2 = new IloRange[it.numC2];

            for (int i : it.C2) {
                if (IntV[i] == 2) {
                    IPcst2[i-it.numC1-1] = IP.addRange(-Double.MAX_VALUE,1.0);
                }
                else if (IntV[i] == 1) {
                    IPcst2[i-it.numC1-1] = IP.addRange(1,1);
                }
                else if (IntV[i] == 0) {
                    IPcst2[i-it.numC1-1] = IP.addRange(-Double.MAX_VALUE,0);
                }
            }

            IloRange IPcst3 = IP.addRange(0.0,it.K);

            ColGen.IloNumVarArray IP_s = new ColGen.IloNumVarArray();

            for (scheme s : schemes) {
                IloColumn cols = IP.column(IPobj, s.getProfit(it));
                int cstInd;
                for (Integer i : s.customers) {
                    if (i >= 1 && i < it.C2[0]) {
                        cstInd = i-1;
                        cols = cols.and(IP.column(IPcst1[cstInd],1));
                    }
                    else {
                        cstInd = i-it.numC1-1;
                        cols = cols.and(IP.column(IPcst2[cstInd],1));
                    }
                }
                cols = cols.and(IP.column(IPcst3,1.0));
                IP_s.add(IP.boolVar(cols));
            }

            IP.setParam(IloCplex.IntParam.RootAlgorithm, IloCplex.Algorithm.Primal);
            IP.setOut(null);

            if (!IP.solve()) {
                System.out.println("er: IMP infeasible!");
                return 1E10;
            }

            for (int i = 0; i < IP_s.getSize(); i++) {
                schemes.get(i).setV(IP.getValue(IP_s.getVar(i)));
                if (Math.abs(schemes.get(i).v - 1) < 1E-7) {
                    schemes.get(i).setV(1);
                }
            }
            IP_obj = IP.getObjValue();
            IP.end();
            return IP_obj;
        } catch (IloException e) {
            System.err.println("Concert exception caught '" + e + "' caught");
        }
        return 1E10;
    }
}
