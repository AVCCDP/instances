import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;

public class InitialSolution {
    ArrayList<scheme> initScheme;

    public InitialSolution(Instance it) {
        this.initScheme = new ArrayList<scheme>();
        for (int i = 1; i < it.numN-1; i++) {
            scheme initscheme = new scheme();
            route initR = new route();
            initR.addcity(0,it);
            initR.addcity(i,it);
            initR.addcity(it.numN-1,it);
            initscheme.setR(initR);
            this.initScheme.add(initscheme);
        }
    }

    public InitialSolution(Instance it){
        this.initScheme = new ArrayList<scheme>();
        ArrayList<Integer> unUseV = new ArrayList<>();
        for (int i = 1; i < it.numN-1; i++) {
            unUseV.add(i);
        }
        int k = 0;
        while(k < it.K && !unUseV.isEmpty()){
            scheme initS = new scheme();
            route initR = new route();
            int lastV = 0;
            double lastQ = 0.0;
            double lastT = 0.0;
            initR.addcity(0,it);
            ArrayList<Integer> KunUseV = new ArrayList<>();
            for (Integer i : unUseV) {
                KunUseV.add(i);
            }
            while ((lastT + it.T[lastV][it.numN-1] <= it.Tmax) && lastQ<=it.W && !KunUseV.isEmpty()) {
                ArrayList<Double> wCoef = new ArrayList<>();
                for (Integer v : KunUseV) {
                    wCoef.add(it.P[v] / (it.T[lastV][v] + it.T[v][it.numN-1] - it.T[lastV][it.numN-1]));
                }
                double maxW = Collections.max(wCoef);
                int maxVid = wCoef.indexOf(maxW);
                int maxV = KunUseV.get(maxVid);
                if ((lastT + it.T[lastV][maxV] + it.T[maxV][it.numN-1] <= it.Tmax) && (lastQ+it.q[maxV] < it.W) && (it.P[maxV] - it.d[maxV]*(lastT + it.T[lastV][maxV]) > 0)) {
                    initR.addcity(maxV,it);
                    lastT = lastT + it.T[lastV][maxV];
                    lastV = maxV;
                    lastQ = lastQ+it.q[maxV];
                    KunUseV.remove(maxVid);
                    unUseV.remove((Integer) maxV);
                } else {
                    KunUseV.remove(maxVid);
                }
            }
            initR.addcity(it.numN-1,it);
            initS.setR(initR);
            initR.setRQ(lastQ);
            initS.setQ(lastQ);
            this.initScheme.add(initS);
            k += 1;
        }

        ArrayList<Integer> movableV = new ArrayList<>();
        for (Integer v : unUseV) {
            ArrayList<Double> allMdW = new ArrayList<>();
            ArrayList<Double> mdW = new ArrayList<>();
            for (int i = 1; i < it.numN-1; i++) {
                allMdW.add(it.md[v][i]);
                mdW.add(it.md[v][i]);
            }
            while (!mdW.isEmpty()) {
                double maxMd = Collections.max(mdW);
                int maxM = allMdW.indexOf(maxMd) + 1;
                for (scheme s : this.initScheme) {
                    int maxMid = s.r.path.indexOf(maxM);
                    if (maxMid > 0 && it.T[v][maxM] <= s.r.pathAT.get(maxMid) && (s.Q + it.q[v]<=it.W) && (it.P[v]*maxMd - it.d[v]*s.r.pathAT.get(maxMid) > 0)) {
                        boolean flag = true;
                        for (move m : s.mList) {
                            if (m.group.get(0) == maxM) {
                                m.addcity(v,it);
                                s.updateCust();
                                s.setQ(s.Q + it.q[v]);
                                m.updateAT(s.r,it,s);
                                movableV.add(v);
                                flag = false;
                            }
                        }
                        if (flag) {
                            move newM = new move();
                            newM.addcity(maxM,it);
                            newM.addcity(v,it);
                            s.setM(newM);
                            s.setQ(s.Q + it.q[v]);
                            newM.updateAT(s.r,it,s);
                            movableV.add(v);
                        }
                        mdW.clear();
                    }
                }
                mdW.remove((Double) maxMd);
            }
        }
        unUseV.removeAll(movableV);

        if (!unUseV.isEmpty()) {
            for (Integer uc : unUseV) {
                if (uc < it.C2[0]) {
                    if (k < it.K) {
                        scheme initscheme = new scheme();
                        route initR = new route();
                        initR.addcity(0,it);
                        initR.addcity(uc,it);
                        initR.addcity(it.numN-1,it);
                        initR.setRQ(it.q[uc]);
                        initscheme.setR(initR);
                        this.initScheme.add(initscheme);
                        k += 1;
                    }
                    else {
                        int sInd = -1;
                        double minNewT = 10000;
                        for (int i = 0; i < this.initScheme.size(); i++) {
                            int lastC = this.initScheme.get(i).r.path.get(this.initScheme.get(i).r.path.size()-2);
                            double lastT = this.initScheme.get(i).r.pathAT.get(this.initScheme.get(i).r.pathAT.size()-2);
                            double newT = lastT + it.T[uc][lastC];
                            if (newT < minNewT) {
                                minNewT = newT;
                                sInd = i;
                            }
                        }

                        route newR = this.initScheme.get(sInd).r;
                        newR.path.add(newR.path.size()-1, uc);
                        newR.pathAT.add(newR.pathAT.size()-1, minNewT);
                        newR.pathLT.add(newR.pathLT.size()-1, minNewT);
                        newR.setRQ(newR.rQ+it.q[uc]);
                        this.initScheme.get(sInd).setQ(newR.rQ);
                        this.initScheme.get(sInd).setR(newR);
                        this.initScheme.get(sInd).customers = new ArrayList<>(new LinkedHashSet<>(this.initScheme.get(sInd).customers));
                    }
                }
            }
        }
    }


    public ArrayList<scheme> getScheme() {
        return this.initScheme;
    }
}
