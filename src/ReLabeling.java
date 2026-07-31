import java.util.*;

public class ReLabeling {
    Instance it;
    ArrayList<rLabel> rLabels;
    ArrayList<route> negRoutes;
    ArrayList<route> posRoutes;
    mLabel[][] mLabels;
    
    class rLabel {
        public int curN;
        public int firstN;
        public int preLid;
        public double rAt;
        public double rQ;
        public double Cr;
        public double lbC;
        public ArrayList<Integer> Path;
        public ArrayList<Double> Time;
        public boolean[] Visited;
        public boolean dominated;
        rLabel(int i, int fi,int preid, double at, double q, double c, double lbc, ArrayList<Integer> p,ArrayList<Double> t, boolean[] v, boolean d){
            curN = i;
            firstN = fi;
            preLid = preid;
            rAt = at;
            rQ = q;
            Cr = c;
            lbC = lbc;
            Path = p;
            Time = t;
            Visited = v;
            dominated = d;
        }
    }
    
    class rLabelComparator implements Comparator<Integer> {
        public int compare (Integer a, Integer b) {
            rLabel A = rLabels.get(a);
            rLabel B = rLabels.get(b);
            if (A.Cr - B.Cr < -1e-7) {
                return -1;
            } else if (A.Cr - B.Cr > 1e-7) {
                return 1;
            } else {
                if (A.curN == B.curN) {
                    if (A.rAt - B.rAt < -1e-7) {
                        return -1;
                    } else if (A.rAt - B.rAt > 1e-7) {
                        return 1;
                    } else {
                        int i = 0;
                        while ( i < it.numN) {
                            if (A.Visited[i] != B.Visited[i]) {
                                if (A.Visited[i]) {
                                    return -1;
                                }else {
                                    return 1;
                                }
                            }
                            i++;
                        }
                        return -1;
                    }
                } else if (A.curN > B.curN) {
                    return 1;
                } else {
                    return -1;
                }
            }
        }
    }
    
    class mLabel{
        public int rendezvous;
        public int moving;
        public double mAt;
        public double mQ;
        public double Cm;
        mLabel(int rv, int m, double at, double q, double c){
            rendezvous = rv;
            moving = m;
            mAt = at;
            mQ = q;
            Cm = c;
        }
    }

    
    public void GenMove(Instance it) {
        int i,j;
        double cm;
        mLabel curLabel;
        mLabels = new mLabel[it.numN][it.numN];
        for (i = 1; i < it.numN-1; i++) {
            for (j = 1; j < it.numN-1; j++) {
                if ((it.feaMove[i][j] == 1) && (it.Tmb[j][i] < it.M - 1e-6)) {
                    cm = it.P_m[j] - it.P[j]*it.md[j][i] + it.Tmb[j][i]*it.d[j];
                    curLabel = new mLabel(i,j,it.Tmb[j][i],it.q[j],cm);
                    mLabels[i][j] = curLabel;
                }
            }
        }
    }

    
    public void GenRouteV3(Instance it, int endNum) {
        rLabel curL;
        int i, j, newid, negRNum,posRNum,maxRNum;
        Integer curLid;
        double newQ;
        double newlbC;
        double newCr;
        double newAt,feaAt;
        ArrayList<Integer> newPath;
        ArrayList<Double> newTime;
        this.it = it;
        this.negRoutes = new ArrayList<route>();
        this.posRoutes = new ArrayList<route>();

        TreeSet<Integer> U = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> negP = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> posP = new TreeSet<Integer>(new rLabelComparator());

        rLabels = new ArrayList<rLabel>(it.numN*100);
        ArrayList<Double> t = new ArrayList<>();
        ArrayList<Integer> p = new ArrayList<>();
        t.add(0.0);
        p.add(0);
        boolean[] visit = new boolean[it.numN];
        visit[0] = true;
        rLabels.add(new rLabel(0,0,-1,0.0,0.0,it.ps,0.0,p,t,visit,false));
        U.add(0);
        ArrayList<Integer>[] allLx = new ArrayList[it.numN];
        for (i = 0; i < it.numN; i++) {
            allLx[i] = new ArrayList<Integer>();
        }
        allLx[0].add(0);
        negRNum = 0;
        posRNum = 0;
        maxRNum = 2 * endNum;
        
        while ((U.size() > 0) && (negRNum < maxRNum) && (posRNum < 400*endNum*endNum)) {
            curLid = U.pollFirst();
            curL = rLabels.get(curLid);
            ArrayList<Integer> Out = new ArrayList<Integer>();
            if (curL.curN != it.numN-1){
                for (i = 0; i < allLx[curL.curN].size(); i++) {
                    int checkid = allLx[curL.curN].get(i);
                    if (checkid != curLid) {
                        rLabel checkL = rLabels.get(checkid);
                        if((curL.rQ <= checkL.rQ) && (curL.lbC <= checkL.lbC) && (curL.rAt <= checkL.rAt)) {
                            rLabels.get(checkid).dominated = true;
                            Out.add(checkid);
                            U.remove(checkid);
                        }
                    }
                }
                for (Integer outid : Out) {
                    allLx[curL.curN].remove(outid);
                }
            }

            if (!curL.dominated) {
                if (curL.curN == it.numN-1) {
                    if (curL.Cr < -1e-7) {
                        negP.add(curLid);
                        negRNum++;
                    } else {
                        posP.add(curLid);
                        posRNum++;
                    }
                }
                else {
                    for (int c : it.C) {
                        newAt = curL.rAt +it.Trb[curL.curN][c];
                        feaAt = newAt + it.Trb[c][it.numN-1];
                        if ((!curL.Visited[c]) && (feaAt <= it.Tmax) && (curL.rQ+it.q[c] <= it.W) && (it.Trb[curL.curN][c] < it.M - 1e-6)) {
                            newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                            newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                            newQ = curL.rQ + it.q[c];
                            newCr = curL.Cr + it.P_r[c] + newAt*it.d[c];
                            newlbC = newCr;
                            for (int uc : it.C) {
                                if (!curL.Visited[uc] && c != uc && it.Trb[c][uc] < it.M - 1e-6) {
                                    double rMinC = 0.0;
                                    if (newAt+it.Trb[c][uc]+it.Trb[uc][it.numN-1] < it.Tmax) {
                                        rMinC = it.P_r[uc] + (newAt+it.Trb[c][uc])*it.d[uc];
                                    }
                                    rMinC = Math.min(rMinC,0);

                                    double minAT = 0;
                                    if (curL.curN == 0) {
                                        minAT = Math.max(it.Trb[0][c], it.minMoveT[uc]);
                                    } else {
                                        minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                    }
                                    double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                    mMinC = Math.min(mMinC,0);
                                    newlbC += Math.min(rMinC,mMinC);
                                }
                            }
                            
                            if ( newlbC < 0 && (newCr <= curL.Cr || newCr < 0)) {
                                newid = rLabels.size();
                                boolean[] newV = new boolean[it.numN];
                                System.arraycopy(curL.Visited,0,newV,0,it.numN);
                                newV[c] = true;
                                newPath.add(c);
                                newTime.add(newAt);
                                if (curL.curN == 0) {
                                    curL.firstN = c;
                                }
                                for (j = 1; j < it.numN-1; j++) {
                                    if (!newV[j]) {
                                        if (newAt + it.Trb[c][j] + it.Trb[j][it.numN-1] > it.Tmax) {
                                            newV[j] = true;
                                        }
                                    }
                                }
                                rLabels.add(new rLabel(c, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));
                                if (!U.add(newid)) {
                                    System.out.println("不能添加");
                                }
                                allLx[c].add(newid);
                            }
                        }
                    }

                    if (curL.curN != 0 && it.Trb[curL.curN][it.numN-1] < it.M - 1e-6) {
                        newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                        newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                        newAt = curL.rAt + it.Trb[curL.curN][it.numN-1];
                        newQ = curL.rQ;
                        newCr = curL.Cr;
                        newlbC = newCr;
                        for (int uc : it.C) {
                            if (!curL.Visited[uc]) {
                                double minAT = 0;
                                minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                newlbC += Math.min(mMinC,0);
                            }
                        }

                        if (newAt <= it.Tmax && newlbC < 0) {
                            newPath.add(it.numN-1);
                            newTime.add(newAt);
                            newid = rLabels.size();
                            boolean[] newV = new boolean[it.numN];
                            System.arraycopy(curL.Visited,0,newV,0,it.numN);
                            newV[it.numN-1] = true;
                            rLabels.add(new rLabel(it.numN-1, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));
                            if (!U.add(newid)) {
                                System.out.println("不能添加");
                            }
                            allLx[it.numN-1].add(newid);
                        }
                    }
                }
            }
        }
        Integer minLID;

        while( (minLID = negP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            for (Integer city : minL.Path) {
                newRoute.addcity(city, it);
            }
            newRoute.setInspectR(minL.Cr);
            newRoute.setRQ(minL.rQ);
            this.negRoutes.add(newRoute);
        }

        while ((minLID = posP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            double checkCr = minL.Cr;
            for (int uc : it.C) {
                if (!minL.Path.contains(uc)) {
                    double minC = 0;
                    for (Integer meet : minL.Path) {
                        if (mLabels[meet][uc] != null && mLabels[meet][uc].Cm < minC) {
                            minC = mLabels[meet][uc].Cm;
                        }
                    }
                    checkCr += minC;
                }
            }
            
            if (checkCr < 0) {
                for (Integer city : minL.Path) {
                    newRoute.addcity(city, it);
                }
                newRoute.setInspectR(minL.Cr);
                newRoute.setRQ(minL.rQ);
                this.posRoutes.add(newRoute);
            }
        }
        it.total_Label += rLabels.size();
        it.countPath += (this.negRoutes.size()+this.posRoutes.size());
        rLabels = null;
    }

    
    public void GenRouteV2(Instance it, int endNum) {
        rLabel curL;
        int i, j, newid, negRNum,posRNum,maxRNum;
        Integer curLid;
        double newQ;
        double newlbC;
        double newCr;
        double newAt,feaAt;
        ArrayList<Integer> newPath;
        ArrayList<Double> newTime;
        this.it = it;
        this.negRoutes = new ArrayList<route>();
        this.posRoutes = new ArrayList<route>();

        TreeSet<Integer> U = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> negP = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> posP = new TreeSet<Integer>(new rLabelComparator());

        rLabels = new ArrayList<rLabel>(it.numN*100);
        ArrayList<Double> t = new ArrayList<>();
        ArrayList<Integer> p = new ArrayList<>();
        t.add(0.0);
        p.add(0);
        boolean[] visit = new boolean[it.numN];
        visit[0] = true;
        rLabels.add(new rLabel(0,0,-1,0.0,0.0,it.ps,0.0,p,t,visit,false));
        U.add(0);

        ArrayList<Integer>[] allLx = new ArrayList[it.numN];
        for (i = 0; i < it.numN; i++) {
            allLx[i] = new ArrayList<Integer>();
        }
        allLx[0].add(0);
        negRNum = 0;
        posRNum = 0;
        maxRNum = 2 * endNum;

        while ((U.size() > 0) && (negRNum < maxRNum)) {
            curLid = U.pollFirst();
            curL = rLabels.get(curLid);
            ArrayList<Integer> Out = new ArrayList<Integer>();
            if (curL.curN != it.numN-1){
                for (i = 0; i < allLx[curL.curN].size(); i++) {
                    int checkid = allLx[curL.curN].get(i);
                    if (checkid != curLid) {
                        rLabel checkL = rLabels.get(checkid);
                        boolean checkDom = true;
                        for (j = 1; j < it.numN-1; j++) {
                            if (curL.Visited[j] && !checkL.Visited[j]) {
                                checkDom = false;
                                break;
                            }
                        }
                        if((curL.rQ <= checkL.rQ) && (curL.lbC <= checkL.lbC) && (curL.rAt <= checkL.rAt)) {
                            rLabels.get(checkid).dominated = true;
                            Out.add(checkid);
                            U.remove(checkid);
                        }
                    }
                }
                for (Integer outid : Out) {
                    allLx[curL.curN].remove(outid);
                }
            }

            if (!curL.dominated) {
                if (curL.curN == it.numN-1) {
                    if (curL.Cr < -1e-7) {
                        negP.add(curLid);
                        negRNum++;
                    } else {
                        posP.add(curLid);
                        posRNum++;
                    }
                }
                else {
                    for (int c : it.C) {
                        if ((!curL.Visited[c]) && (curL.rQ+it.q[c] <= it.W) && (it.Trb[curL.curN][c] < it.M - 1e-6)) {
                            newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                            newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                            newAt = curL.rAt +it.Trb[curL.curN][c];
                            feaAt = newAt + it.Trb[c][it.numN-1];
                            newQ = curL.rQ + it.q[c];
                            newCr = curL.Cr + it.P_r[c] + newAt*it.d[c];
                            newlbC = newCr;
                            for (int uc : it.C) {
                                if (!curL.Path.contains(uc) && c != uc && it.Trb[c][uc] < it.M - 1e-6) {
                                    double rMinC = 0.0;
                                    if (newAt+it.Trb[c][uc]+it.Trb[uc][it.numN-1] < it.Tmax) {
                                        rMinC = it.P_r[uc] + (newAt+it.Trb[c][uc])*it.d[uc];
                                    }
                                    rMinC = Math.min(rMinC,0);

                                    double minAT = 0;
                                    if (curL.curN == 0) {
                                        minAT = Math.max(it.Trb[0][c], it.minMoveT[uc]);
                                    } else {
                                        minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                    }
                                    double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                    mMinC = Math.min(mMinC,0);
                                    newlbC += Math.min(rMinC,mMinC);
                                }
                            }

                            if (feaAt <= it.Tmax && newQ <= it.W && newlbC < 0) {
                                newid = rLabels.size();
                                boolean[] newV = new boolean[it.numN];
                                System.arraycopy(curL.Visited,0,newV,0,it.numN);
                                newV[c] = true;
                                newPath.add(c);
                                newTime.add(newAt);
                                if (curL.curN == 0) {
                                    curL.firstN = c;
                                }
                                
                                for (j = 1; j < it.numN-1; j++) {
                                    if (!newV[j]) {
                                        if (newAt + it.Trb[c][j] + it.Trb[j][it.numN-1] > it.Tmax) {
                                            newV[j] = true;
                                        }
                                    }
                                }
                                
                                rLabels.add(new rLabel(c, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));

                                if (!U.add(newid)) {
                                    System.out.println("不能添加");
                                }
                                allLx[c].add(newid);
                            }
                        }
                    }

                    if (curL.curN != 0 && it.Trb[curL.curN][it.numN-1] < it.M - 1e-6) {
                        newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                        newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                        newAt = curL.rAt + it.Trb[curL.curN][it.numN-1];
                        newQ = curL.rQ;
                        newCr = curL.Cr;
                        newlbC = newCr;
                        for (int uc : it.C) {
                            if (!curL.Path.contains(uc)) {
                                double minAT = 0;
                                minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                newlbC += Math.min(mMinC,0);
                            }
                        }

                        if (newAt <= it.Tmax && newlbC < 0) {
                            newid = rLabels.size();
                            boolean[] newV = new boolean[it.numN];
                            System.arraycopy(curL.Visited,0,newV,0,it.numN);
                            newV[it.numN-1] = true;
                            newPath.add(it.numN-1);
                            newTime.add(newAt);
                            rLabels.add(new rLabel(it.numN-1, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));
                            if (!U.add(newid)) {
                                System.out.println("不能添加");
                            }
                            allLx[it.numN-1].add(newid);
                        }
                    }
                }
            }
        }

        Integer minLID;

        while( (minLID = negP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            for (Integer city : minL.Path) {
                newRoute.addcity(city, it);
            }
            newRoute.setInspectR(minL.Cr);
            newRoute.setRQ(minL.rQ);
            this.negRoutes.add(newRoute);
        }

        while ((minLID = posP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            double checkCr = minL.Cr;
            for (int uc : it.C) {
                if (!minL.Path.contains(uc)) {
                    double minC = 0;
                    for (Integer meet : minL.Path) {
                        if (mLabels[meet][uc] != null && mLabels[meet][uc].Cm < minC) {
                            minC = mLabels[meet][uc].Cm;
                        }
                    }
                    checkCr += minC;
                }
            }
            
            if (checkCr < 0) {
                for (Integer city : minL.Path) {
                    newRoute.addcity(city, it);
                }
                newRoute.setInspectR(minL.Cr);
                newRoute.setRQ(minL.rQ);
                this.posRoutes.add(newRoute);
            }
        }
        
        it.total_Label += rLabels.size();
        rLabels = null;
    }

    
    public void GenRouteV1(Instance it, int endNum) {
        rLabel curL;
        int i, j, newid, negRNum,posRNum,maxRNum;
        Integer curLid;
        double newQ;
        double newlbC;
        double newCr;
        double newAt,feaAt;
        ArrayList<Integer> newPath;
        ArrayList<Double> newTime;
        this.it = it;
        this.negRoutes = new ArrayList<route>();
        this.posRoutes = new ArrayList<route>();

        TreeSet<Integer> U = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> negP = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> posP = new TreeSet<Integer>(new rLabelComparator());

        rLabels = new ArrayList<rLabel>(it.numN*100);
        ArrayList<Double> t = new ArrayList<>();
        ArrayList<Integer> p = new ArrayList<>();
        t.add(0.0);
        p.add(0);
        boolean[] visit = new boolean[it.numN];
        visit[0] = true;
        rLabels.add(new rLabel(0,0,-1,0.0,0.0,it.ps,0.0,p,t,visit,false));
        U.add(0);
        ArrayList<Integer>[] allLx = new ArrayList[it.numN];
        for (i = 0; i < it.numN; i++) {
            allLx[i] = new ArrayList<Integer>();
        }
        allLx[0].add(0);
        negRNum = 0;
        posRNum = 0;
        maxRNum = 2 * endNum;
        
        while ((U.size() > 0) && (negRNum < maxRNum)) {
            if (System.currentTimeMillis() - it.startT >= 7200*1000) {
                System.out.println("Time limit, end");
                it.overTime = true;
                break;
            }
            curLid = U.pollFirst();
            curL = rLabels.get(curLid);
            ArrayList<Integer> Out = new ArrayList<Integer>();
            
            if (!curL.dominated) {
                if (curL.curN == it.numN-1) {
                    if (curL.Cr < -1e-7) {
                        negP.add(curLid);
                        negRNum++;
                    } else {
                        posP.add(curLid);
                        posRNum++;
                    }
                }
                else {
                    for (int c : it.C) {
                        newAt = curL.rAt +it.Trb[curL.curN][c];
                        feaAt = newAt + it.Trb[c][it.numN-1];
                        if ((!curL.Visited[c]) &&(feaAt<=it.Tmax) && (curL.rQ+it.q[c] <= it.W) && (it.Trb[curL.curN][c] < it.M - 1e-6)) {
                            newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                            newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                            newQ = curL.rQ + it.q[c];
                            newCr = curL.Cr + it.P_r[c] + newAt*it.d[c];
                            newlbC = newCr;
                            for (int uc : it.C) {
                                if (!curL.Visited[uc] && c != uc && it.Trb[c][uc] < it.M - 1e-6) {
                                    double rMinC = 0.0;
                                    if (newAt+it.Trb[c][uc]+it.Trb[uc][it.numN-1] < it.Tmax) {
                                        rMinC = it.P_r[uc] + (newAt+it.Trb[c][uc])*it.d[uc];
                                    }
                                    rMinC = Math.min(rMinC,0);

                                    double minAT = 0;
                                    if (curL.curN == 0) {
                                        minAT = Math.max(it.Trb[0][c], it.minMoveT[uc]);
                                    } else {
                                        minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                    }
                                    double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                    mMinC = Math.min(mMinC,0);
                                    newlbC += Math.min(rMinC,mMinC);
                                }
                            }


                            if ( newlbC < 0) {
                                newid = rLabels.size();
                                boolean[] newV = new boolean[it.numN];
                                System.arraycopy(curL.Visited,0,newV,0,it.numN);
                                newV[c] = true;
                                newPath.add(c);
                                newTime.add(newAt);
                                if (curL.curN == 0) {
                                    curL.firstN = c;
                                }

                                rLabels.add(new rLabel(c, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime,newV, false));
                                if (!U.add(newid)) {
                                    System.out.println("不能添加");
                                }
                                allLx[c].add(newid);
                            }
                        }
                    }

                    if (curL.curN != 0 && it.Trb[curL.curN][it.numN-1] < it.M - 1e-6) {
                        newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                        newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                        newAt = curL.rAt + it.Trb[curL.curN][it.numN-1];
                        newQ = curL.rQ;
                        newCr = curL.Cr;
                        newlbC = newCr;
                        for (int uc : it.C) {
                            if (!curL.Visited[uc]) {
                                double minAT = 0;
                                minAT = Math.max(it.Trb[0][curL.firstN], it.minMoveT[uc]);
                                double mMinC = it.P_m[uc] - it.P[uc]*it.mdMax[uc] + minAT*it.d[uc];
                                newlbC += Math.min(mMinC,0);
                            }
                        }

                        if (newAt <= it.Tmax && newlbC < 0) {
                            newid = rLabels.size();
                            boolean[] newV = new boolean[it.numN];
                            System.arraycopy(curL.Visited,0,newV,0,it.numN);
                            newV[it.numN-1] = true;
                            newPath.add(it.numN-1);
                            newTime.add(newAt);
                            rLabels.add(new rLabel(it.numN-1, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime,newV, false));
                            if (!U.add(newid)) {
                                System.out.println("不能添加");
                            }
                            allLx[it.numN-1].add(newid);
                        }
                    }
                }
            }
        }


        Integer minLID;

        while( (minLID = negP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            for (Integer city : minL.Path) {
                newRoute.addcity(city, it);
            }
            newRoute.setInspectR(minL.Cr);
            newRoute.setRQ(minL.rQ);
            this.negRoutes.add(newRoute);
        }

        while ((minLID = posP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            double checkCr = minL.Cr;
            for (int uc : it.C) {
                if (!minL.Path.contains(uc)) {
                    double minC = 0;
                    for (Integer meet : minL.Path) {
                        if (mLabels[meet][uc] != null && mLabels[meet][uc].Cm < minC) {
                            minC = mLabels[meet][uc].Cm;
                        }
                    }
                    checkCr += minC;
                }
            }

            if (checkCr < 0) {
                for (Integer city : minL.Path) {
                    newRoute.addcity(city, it);
                }
                newRoute.setInspectR(minL.Cr);
                newRoute.setRQ(minL.rQ);
                this.posRoutes.add(newRoute);
            }
        }
        it.total_Label += rLabels.size();
        it.countPath += (this.negRoutes.size()+this.posRoutes.size());
        rLabels = null;
    }

    
    public void GenRouteV0(Instance it, int endNum) {
        rLabel curL;
        int i, j, newid, negRNum,posRNum,maxRNum;
        Integer curLid;
        double newQ;
        double newlbC;
        double newCr;
        double newAt,feaAt;
        ArrayList<Integer> newPath;
        ArrayList<Double> newTime;
        this.it = it;
        this.negRoutes = new ArrayList<route>();
        this.posRoutes = new ArrayList<route>();

        TreeSet<Integer> U = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> negP = new TreeSet<Integer>(new rLabelComparator());

        TreeSet<Integer> posP = new TreeSet<Integer>(new rLabelComparator());

        rLabels = new ArrayList<rLabel>(it.numN*100);
        ArrayList<Double> t = new ArrayList<>();
        ArrayList<Integer> p = new ArrayList<>();
        t.add(0.0);
        p.add(0);
        boolean[] visit = new boolean[it.numN];
        visit[0] = true;
        rLabels.add(new rLabel(0,0,-1,0.0,0.0,it.ps,0.0,p,t,visit,false));
        U.add(0);
        ArrayList<Integer>[] allLx = new ArrayList[it.numN];
        for (i = 0; i < it.numN; i++) {
            allLx[i] = new ArrayList<Integer>();
        }
        allLx[0].add(0);
        negRNum = 0;
        posRNum = 0;
        maxRNum = 2 * endNum;

        while ((U.size() > 0) && (negRNum < maxRNum)) {
            curLid = U.pollFirst();
            curL = rLabels.get(curLid);
            if (!curL.dominated) {
                if (curL.curN == it.numN-1) {
                    if (curL.Cr < -1e-7) {
                        negP.add(curLid);
                        negRNum++;
                    } else {
                        posP.add(curLid);
                        posRNum++;
                    }
                }
                else {
                    for (int c : it.C) {
                        if ((!curL.Visited[c]) && (curL.rQ+it.q[c] <= it.W) && (it.Trb[curL.curN][c] < it.M - 1e-6)) {
                            newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                            newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                            newAt = curL.rAt +it.Trb[curL.curN][c];
                            feaAt = newAt + it.Trb[c][it.numN-1];
                            newQ = curL.rQ + it.q[c];
                            newCr = curL.Cr + it.P_r[c] + newAt*it.d[c];
                            newlbC = newCr;

                            if (feaAt <= it.Tmax && newQ <= it.W ) {
                                newid = rLabels.size();
                                boolean[] newV = new boolean[it.numN];
                                System.arraycopy(curL.Visited,0,newV,0,it.numN);
                                newV[c] = true;
                                newPath.add(c);
                                newTime.add(newAt);
                                if (curL.curN == 0) {
                                    curL.firstN = c;
                                }
                                for (j = 1; j < it.numN-1; j++) {
                                    if (!newV[j]) {
                                        if (newAt + it.Trb[c][j] + it.Trb[j][it.numN-1] > it.Tmax) {
                                            newV[j] = true;
                                        }
                                    }
                                }

                                rLabels.add(new rLabel(c, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));
                                if (!U.add(newid)) {
                                    System.out.println("不能添加");
                                }
                                allLx[c].add(newid);
                            }
                        }
                    }

                    if (curL.curN != 0 && it.Trb[curL.curN][it.numN-1] < it.M - 1e-6) {
                        newPath = new ArrayList<>(curL.Path.subList(0, curL.Path.size()));
                        newTime = new ArrayList<>(curL.Time.subList(0, curL.Time.size()));
                        newAt = curL.rAt + it.Trb[curL.curN][it.numN-1];
                        newQ = curL.rQ;
                        newCr = curL.Cr;
                        newlbC = newCr;
                        if (newAt <= it.Tmax ) {
                            newid = rLabels.size();
                            boolean[] newV = new boolean[it.numN];
                            System.arraycopy(curL.Visited,0,newV,0,it.numN);
                            newV[it.numN-1] = true;
                            newPath.add(it.numN-1);
                            newTime.add(newAt);
                            rLabels.add(new rLabel(it.numN-1, curL.firstN, curLid, newAt, newQ, newCr, newlbC, newPath, newTime, newV, false));
                            if (!U.add(newid)) {
                                System.out.println("不能添加");
                            }
                            allLx[it.numN-1].add(newid);
                        }
                    }
                }
            }
        }


        int MB = 1024 * 1024;
        Runtime runtime = Runtime.getRuntime();
        System.out.println("找完标签后内存：");
        System.out.println("Java Memory=> Total:" + (runtime.totalMemory() / MB)
                + " Max:" + (runtime.maxMemory() / MB) + " Used:"
                + ((runtime.totalMemory() - runtime.freeMemory()) / MB) + " Free: "
                + runtime.freeMemory() / MB);

        Integer minLID;

        while( (minLID = negP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            for (Integer city : minL.Path) {
                newRoute.addcity(city, it);
            }
            newRoute.setInspectR(minL.Cr);
            newRoute.setRQ(minL.rQ);
            this.negRoutes.add(newRoute);
        }

        while ((minLID = posP.pollFirst()) != null) {
            rLabel minL = rLabels.get(minLID);
            route newRoute = new route();
            double checkCr = minL.Cr;
            for (int uc : it.C) {
                if (!minL.Path.contains(uc)) {
                    double minC = 0;
                    for (Integer meet : minL.Path) {
                        if (mLabels[meet][uc] != null && mLabels[meet][uc].Cm < minC) {
                            minC = mLabels[meet][uc].Cm;
                        }
                    }
                    checkCr += minC;
                }
            }

            if (checkCr < 0) {
                for (Integer city : minL.Path) {
                    newRoute.addcity(city, it);
                }
                newRoute.setInspectR(minL.Cr);
                newRoute.setRQ(minL.rQ);
                this.posRoutes.add(newRoute);
            }
        }

        it.total_Label += rLabels.size();
        it.countPath += (this.negRoutes.size()+this.posRoutes.size());
        rLabels = null;
    }
    
    public void GenScheme(ArrayList<scheme> newSchemes) throws CloneNotSupportedException {
        int i,j;
        ArrayList<Integer> SelectM = new ArrayList<Integer>(it.numC);
        ArrayList<Integer> SelectC = new ArrayList<Integer>(it.numC);

        if (!this.negRoutes.isEmpty()) {
            for (route r : negRoutes) {
                route r0 = r.clone();
                scheme newS = new scheme();
                double newQ = r.rQ;
                boolean flag = false;
                mLabels = new mLabel[it.numN][it.numN];
                mLabel curLabel;
                double cm;
                double atm;
                int numMove = 0;
                for (Integer meet : r.path) {
                    for (int c : it.C) {
                        if ( !r.path.contains(c) && r.rQ+it.q[c] <= it.W && it.feaMove[meet][c] == 1 && it.Tmb[c][meet] < it.M - 1e-6) {
                            atm = Math.max(it.Tmb[c][meet], r.pathAT.get(r.path.indexOf(meet)));
                            cm = it.P_m[c] - it.P[c]*it.md[c][meet] + atm*it.d[c];
                            if (cm < -1*r.inspectR) {
                                curLabel = new mLabel(meet, c, atm, it.q[c], cm);
                                mLabels[meet][c] = curLabel;
                                numMove += 1;
                            }
                        }
                    }
                }

                while (numMove > 0) {
                    route tempR = r.clone();
                    scheme tempS = newS.clone();
                    double minCm = -r.inspectR;
                    int minM = -1;
                    int minC = -1;
                    for (i = 1; i < it.numN-1; i++) {
                        for (j = 1; j < it.numN-1; j++) {
                            if (mLabels[i][j] != null && mLabels[i][j].Cm < minCm) {
                                minCm = mLabels[i][j].Cm;
                                minM = i;
                                minC = j;
                            }
                        }
                    }

                    int Mid = SelectM.indexOf(minM);

                    if (minC!= -1 && Mid == -1 && !SelectC.contains(minC)) {
                        move m = new move();
                        m.addcity(minM,it);
                        m.addcity(minC,it);
                        m.updateAT(tempR,it,tempS);
                        if (tempR.pathAT.get(tempR.path.size()-1) <= it.Tmax && newQ+it.q[minC] <= it.W) {
                            double inspectSol = 0.0;
                            for (move okm : tempS.mList) {
                                inspectSol += okm.getInspectM(it);
                            }
                            inspectSol += tempR.getInspectR(it) + m.getInspectM(it);
                            if (inspectSol < -1e-6) {
                                newQ += it.q[minC];
                                tempS.setM(m);
                                newS = tempS.clone();
                                newS.setInspectS(inspectSol);
                                r = tempR.clone();
                                SelectM.add(minM);
                                SelectC.add(minC);
                                flag = true;
                            }
                        }
                    }
                    else if (minC!= -1 && Mid != -1 && !SelectC.contains(minC)) {
                        move tempM = newS.mList.get(Mid).clone();
                        tempM.addcity(minC,it);
                        tempM.updateAT(tempR,it,tempS);
                        if (tempR.pathAT.get(tempR.path.size()-1) <= it.Tmax && newQ+it.q[minC] <= it.W){
                            double inspectSol = 0.0;
                            for (int mi = 0; mi < tempS.mList.size(); mi++) {
                                if (mi != Mid) {
                                    inspectSol += tempS.mList.get(mi).getInspectM(it);
                                }
                            }

                            inspectSol += tempR.getInspectR(it) + tempM.getInspectM(it);

                            if (inspectSol < -1e-6) {
                                newQ += it.q[minC];
                                tempS.mList.set(Mid,tempM);
                                newS = tempS.clone();
                                newS.updateCust();
                                newS.setInspectS(inspectSol);
                                r = tempR.clone();
                                SelectC.add(minC);
                                flag = true;
                            }
                        }
                    }

                    if (minM != -1 && minC != -1) {
                        mLabels[minM][minC] = null;
                    }
                    numMove -= 1;
                }
                if (flag) {
                    newS.setR(r);
                    newS.setQ(newQ);

                    newSchemes.add(newS);
                }
                scheme new2S = new scheme();
                new2S.setR(r0);
                new2S.setQ(r0.rQ);
                newSchemes.add(new2S);
                SelectM.clear();
                SelectC.clear();
            }
        }
        else if (!this.posRoutes.isEmpty()) {
            for (route r : posRoutes) {
                scheme newS = new scheme();
                double newQ = r.rQ;
                double ckCost = r.inspectR;
                mLabels = new mLabel[it.numN][it.numN];
                mLabel curLabel;
                double cm;
                double atm;
                int numMove = 0;
                for (Integer meet : r.path) {
                    for (int c : it.C) {
                        if ( meet != 0 && !r.path.contains(c) && r.rQ+it.q[c] <= it.W && it.feaMove[meet][c] == 1 && it.Tmb[c][meet] < it.M - 1e-6) {
                            atm = Math.max(it.Tmb[c][meet], r.pathAT.get(r.path.indexOf(meet)));
                            cm = it.P_m[c] - it.P[c]*it.md[c][meet] + atm*it.d[c];
                            if (cm < -1e-7) {
                                curLabel = new mLabel(meet, c, atm, it.q[c], cm);
                                mLabels[meet][c] = curLabel;
                                numMove += 1;
                            }
                        }
                    }
                }

                while (numMove > 0) {
                    route tempR = r.clone();
                    scheme tempS = newS.clone();
                    double minCm = 0;
                    int minM = -1;
                    int minC = -1;
                    for (i = 1; i < it.numN-1; i++) {
                        for (j = 1; j < it.numN-1; j++) {
                            if (mLabels[i][j] != null && mLabels[i][j].Cm < minCm) {
                                minCm = mLabels[i][j].Cm;
                                minM = i;
                                minC = j;
                            }
                        }
                    }
                    int Mid = SelectM.indexOf(minM);
                    if (minC!= -1 && Mid == -1 && !SelectC.contains(minC)) {
                        move m = new move();
                        m.addcity(minM,it);
                        m.addcity(minC,it);
                        m.updateAT(tempR,it,tempS);
                        if (tempR.pathAT.get(tempR.path.size()-1) <= it.Tmax && newQ+it.q[minC] <= it.W) {
                            double newCKC = 0.0;
                            for (move okm : tempS.mList) {
                                newCKC += okm.getInspectM(it);
                            }
                            newCKC += tempR.getInspectR(it) + m.getInspectM(it);
                            if (newCKC < ckCost || newCKC < -1e-6) {
                                ckCost = newCKC;
                                newQ += it.q[minC];
                                tempS.setM(m);
                                newS = tempS.clone();
                                newS.setInspectS(newCKC);
                                r = tempR.clone();
                                SelectM.add(minM);
                                SelectC.add(minC);
                            }
                        }
                    }
                    else if (minC!= -1 && Mid != -1 && !SelectC.contains(minC)) {
                        move tempM = newS.mList.get(Mid).clone();
                        tempM.addcity(minC,it);
                        tempM.updateAT(tempR,it,tempS);
                        if (tempR.pathAT.get(tempR.path.size()-1) <= it.Tmax && newQ+it.q[minC] <= it.W){
                            double newCKC = 0.0;
                            for (int mi = 0; mi < tempS.mList.size(); mi++) {
                                if (mi != Mid) {
                                    newCKC += tempS.mList.get(mi).getInspectM(it);
                                }
                            }
                            newCKC += tempR.getInspectR(it) + tempM.getInspectM(it);

                            if (newCKC < ckCost || newCKC < -1e-6) {
                                ckCost = newCKC;
                                newQ += it.q[minC];
                                tempS.mList.set(Mid,tempM);
                                newS = tempS.clone();
                                newS.updateCust();
                                newS.setInspectS(newCKC);
                                r = tempR.clone();
                                SelectC.add(minC);
                            }
                        }
                    }
                    mLabels[minM][minC] = null;
                    numMove -= 1;
                }

                if (ckCost < -1e-6) {
                    newS.setR(r);
                    newS.setQ(newQ);
                    newSchemes.add(newS);
                }
                SelectM.clear();
                SelectC.clear();
            }
        }
    }
}
