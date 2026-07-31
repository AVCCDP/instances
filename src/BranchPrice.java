import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BranchPrice {
    double LB = -1E10;
    double UB = 1E10;
    int totalNodes;

    public BranchPrice() {
        this.totalNodes = BBtree.nodeCount;
    }
    
    class BBtree {
        static int nodeCount = 0;
        BBtree father;
        BBtree son0;
        BBtree son1;
        int bArcS;
        int bArcT;
        int branch;
        double nodeUB;
        boolean toplevel;
        int[] IntV;
        boolean rBranch;
        public BBtree() {
            nodeCount++;
        }
        public static int getNodeCount() {
            return nodeCount;
        }
    }

    
    public void arcBranching(Instance it, BBtree bbt, boolean recur) {
        int i;
        if (bbt.father != null) {
            if (bbt.branch == 0 && bbt.rBranch) {
                it.Trb[bbt.bArcS][bbt.bArcT] =it.M;
            }
            else if (bbt.branch == 1 && bbt.rBranch) {
                for (i = 0; i < bbt.bArcT; i++) {
                    it.Trb[bbt.bArcS][i] = it.M;
                }
                for (i++; i < it.numN; i++) {
                    it.Trb[bbt.bArcS][i] = it.M;
                }
                for (i = 0; i < bbt.bArcS; i++) {
                    it.Trb[i][bbt.bArcT] = it.M;
                }
                for (i++; i < it.numN; i++) {
                    it.Trb[i][bbt.bArcT] = it.M;
                }
                it.Trb[bbt.bArcT][bbt.bArcS] = it.M;
                for (i = 0; i < it.numN; i++) {
                    it.Tmb[bbt.bArcS][i] = it.M;
                }
                for (i = 0; i < it.numN; i++) {
                    it.Tmb[bbt.bArcT][i] = it.M;
                }
            }

            else if (bbt.branch == 0 && !bbt.rBranch && bbt.bArcS != -1) {
                it.Tmb[bbt.bArcS][bbt.bArcT] =it.M;
            }

            else if (bbt.branch == 1 && !bbt.rBranch && bbt.bArcS != -1) {
                for (i = 0; i < bbt.bArcT; i++) {
                    it.Tmb[bbt.bArcS][i] = it.M;
                }
                for (i++; i < it.numN; i++) {
                    it.Tmb[bbt.bArcS][i] = it.M;
                }
                for (i = 0; i < it.numN; i++) {
                    it.Tmb[i][bbt.bArcS] = it.M;
                }
                for (i = 0; i < it.numN; i++) {
                    it.Trb[bbt.bArcS][i] = it.M;
                    it.Trb[i][bbt.bArcS] = it.M;
                }
            }
        }

        if (recur && bbt.father.father != null) {
            arcBranching(it, bbt.father,true);
        }
    }

    
    public boolean BBNode(Instance it, ArrayList<scheme> schemes,
                          ArrayList<scheme> bestSchemes, BBtree bbt, int depth)
            throws IOException {
        int i, j, iNode, jNode, selectArc1, selectArc2,branchValue;
        double CGobj,closeness,bestCloseness;
        double[] sCount;
        boolean sInteger,aInteger,isRB;
        try {
            if ((UB - LB) / UB < it.gap) {
                return true;

            if (bbt == null) {
                BBtree newNode = new BBtree();
                newNode.father = null;
                newNode.son0 = null;
                newNode.son1 = null;
                newNode.toplevel = true;
                newNode.bArcS = -1;
                newNode.bArcT = -1;
                newNode.branch = -1;
                newNode.IntV = new int[it.numN];
                Arrays.fill(newNode.IntV,2);
                newNode.rBranch = false;
                bbt = newNode;
            }
            if (bbt.branch < 1) {
                System.out.println("\nforbid: Arc(" + bbt.bArcS + " , " + bbt.bArcT + ")");
            } else {
                System.out.println("\nset: Arc(" + bbt.bArcS + " , " + bbt.bArcT + ")");
            }
            int MB = 1024 * 1024;

            Runtime runtime = Runtime.getRuntime();
            System.out.println("Java Memory=> Total:" + (runtime.totalMemory() / MB)
                    + " Max:" + (runtime.maxMemory() / MB) + " Used:"
                    + ((runtime.totalMemory() - runtime.freeMemory()) / MB) + " Free: "
                    + runtime.freeMemory() / MB);
            ColGen CG = new ColGen();
            CGobj = CG.runColGen(it,schemes,bbt.IntV);


            if (bbt.father == null) {
                UB = CGobj;
            } else if (bbt.father != null && bbt.father.son0 != null && bbt.father.toplevel) {
                UB = Math.min(bbt.nodeUB,bbt.father.nodeUB);
                bbt.toplevel = true;
            }

            if ((CGobj > it.totalP) || (CGobj < 0)) {
                System.out.println("Infeasible | Lower bound: " + LB
                        + " | Upper bound: " + UB + " | Gap: "
                        + ((UB - LB) / UB) + " | BB Depth: "
                        + depth + " | " + "schemes: " + schemes.size());
                return true;
            }

            if (bbt.nodeUB < LB) {
                CG = null;
                System.out.println("Cut | Lower bound: " + LB
                        + " | Upper bound: " + UB + " | Gap: "
                        + ((UB - LB) / UB) + " | BB Depth: "
                        + depth + " | Local CG: " + CGobj + " | " + "schemes: " + schemes.size());
                return true;
            }

            else {
                sInteger = true;
                aInteger = true;
                selectArc1 = -1;
                selectArc2 = -1;
                bestCloseness = -1.0;
                branchValue = 0;
                sCount = new double[it.numN];
                isRB = false;

                for (i = 0; i < it.numN; i++) {
                    java.util.Arrays.fill(it.rArc[i], 0.0);
                    java.util.Arrays.fill(it.mArc[i], 0.0);
                }
                for (scheme s : schemes) {
                    if (s.getV() > 0) {
                        ArrayList<Integer> path = s.r.getpath();
                        iNode = 0;
                        for (i = 1; i < path.size(); i++) {
                            jNode = path.get(i);
                            sCount[jNode] += s.getV();
                            it.rArc[iNode][jNode] += s.getV();
                            iNode = jNode;
                        }
                        for (move m : s.mList) {
                            ArrayList<Integer> group = m.getgroup();
                            iNode = group.get(0);
                            for (i = 1; i < group.size(); i++) {
                                jNode = group.get(i);

                                sCount[jNode] += s.getV();

                                it.mArc[jNode][iNode] += s.getV();
                            }
                        }
                    }
                }
                for (int c : it.C) {
                    if (sCount[c] >= 0.99 ){
                        sCount[c] = 1;
                    }
                }

                int fractInd = 0;
                for (i = 1; i < it.numN-1; i++) {
                    if (Math.floor(sCount[i]) != sCount[i]) {
                        closeness = Math.min(sCount[i],Math.abs(1.0 - sCount[i]));
                        closeness = it.P[i]*closeness;
                        if (closeness > bestCloseness) {
                            fractInd = i;
                            bestCloseness = closeness;
                            sInteger = false;
                        }
                    }
                }
                bbt.IntV[fractInd] = 0;

                if (sInteger) {
                    for (i = 0; i < it.numN; i++) {
                        for (j = 0; j < it.numN; j++) {
                            if ((it.rArc[i][j] > 0) && (it.rArc[i][j] < 0.9999999) || (it.rArc[i][j] >1.00000001)) {
                                if (i > 0 && j <it.numN-1) {
                                    closeness = Math.min(it.rArc[i][j],Math.abs(1.0 - it.rArc[i][j]));
                                    closeness *= (it.P[i]+it.P[j]);
                                    if (closeness > bestCloseness) {
                                        selectArc1 = i;
                                        selectArc2 = j;
                                        bestCloseness = closeness;
                                        branchValue = (Math.abs(1.0 - it.rArc[i][j]) > it.rArc[i][j]) ? 0 : 1;
                                        aInteger = false;
                                        isRB = true;
                                    }
                                }
                            }
                        }
                    }
                    if (aInteger) {
                        for (i = 1;  i< it.numN-1; i++) {
                            for (j = 1;  j< it.numN-1; j++) {
                                if ((it.mArc[i][j] > 0) && ((it.mArc[i][j] < 0.9999999999) || (it.mArc[i][j] > 1.000000001))) {
                                    closeness = Math.min(it.mArc[i][j],Math.abs(1.0 - it.mArc[i][j]));
                                    closeness *= (it.P[i]+it.P[j]);
                                    if (closeness > bestCloseness) {
                                        selectArc1 = i;
                                        selectArc2 = j;
                                        bestCloseness = closeness;
                                        aInteger = false;
                                    }
                                }
                            }
                        }
                    }
                }

                if (sInteger && aInteger) {

                    if (bbt.nodeUB > LB) {
                        LB = bbt.nodeUB;
                        bestSchemes.clear();
                        for (scheme s : schemes) {
                            if (s.getV()>0){
                                bestSchemes.add(s);
                            }
                        }
                        System.out.println("OPT | Lower bound: " + LB
                                + " | Upper bound: " + UB + " | Gap: "
                                + ((UB - LB) / UB) + " | BB Depth: "
                                + depth + " | Local CG: " + CGobj + " | " + "schemes: " + schemes.size());
                    }
                    else {
                        System.out.println("Fea | Lower bound: " + LB
                                + " | Upper bound: " + UB + " | Gap: "
                                + ((UB - LB) / UB) + " | BB Depth: "
                                + depth + " | Local CG: " + CGobj + " | " + "schemes: " + schemes.size());
                    }
                    return true;
                }

                else if (!sInteger) {
                    if (bbt.father == null) {
                        repairMP rMP = new repairMP();
                        rMP.runIP(it,schemes,bbt.IntV);
                        LB = rMP.IP_obj;
                        bestSchemes.clear();
                        for (scheme s : schemes) {
                            if (s.v > 0){
                                bestSchemes.add(s);
                            }
                        }
                        if (it.overTime) {
                            return true;
                        }
                    }
                    System.out.println("Fract | Lower bound: " + LB
                            + " | Upper bound: " + UB + " | Gap: "
                            + ((UB - LB) / UB) + " | BB Depth: "
                            + depth + " | Local CG: " + CGobj + " | " + "schemes: " + schemes.size());
                    System.out.flush();

                    BBtree newNode0 = new BBtree();
                    newNode0.father = bbt;
                    newNode0.son0 = null;
                    newNode0.son1 = null;
                    newNode0.bArcS = selectArc1;
                    newNode0.bArcT = selectArc2;
                    newNode0.branch = branchValue;
                    newNode0.nodeUB = 1E10;
                    newNode0.rBranch = false;
                    newNode0.IntV = new int[it.numN];
                    System.arraycopy(bbt.IntV,0,newNode0.IntV,0,it.numN);
                    ArrayList<scheme> node0Schemes = new ArrayList<scheme>();
                    for (scheme s : schemes) {
                        node0Schemes.add(s);
                    }
                    boolean subTree = BBNode(it,node0Schemes,bestSchemes,newNode0,depth + 1);
                    node0Schemes = null;
                    if (!subTree) {
                        return false;
                    }
                    bbt.son0 = newNode0;

                    BBtree newNode1 = new BBtree();
                    newNode1.father = bbt;
                    newNode1.son0 = null;
                    newNode1.son1 = null;
                    newNode1.bArcS = selectArc1;
                    newNode1.bArcT = selectArc2;
                    newNode1.branch = 1 - branchValue;
                    newNode1.nodeUB =1E10;
                    newNode1.rBranch = false;
                    newNode1.IntV = new int[it.numN];
                    System.arraycopy(bbt.IntV,0,newNode1.IntV,0,it.numN);
                    newNode1.IntV[fractInd] = 1;
                    ArrayList<scheme> node1Schemes = new ArrayList<scheme>();
                    for (scheme s : schemes) {
                        node1Schemes.add(s);
                    }
                    subTree = BBNode(it,node1Schemes,bestSchemes,newNode1,depth + 1);
                    node1Schemes = null;
                    bbt.son1 = newNode1;
                    bbt.nodeUB = Math.max(newNode0.nodeUB, newNode1.nodeUB);
                    return subTree;
                }

                else {
                    if (bbt.father == null) {
                        repairMP rMP = new repairMP();
                        rMP.runIP(it,schemes,bbt.IntV);
                        LB = rMP.IP_obj;
                        bestSchemes.clear();
                        for (scheme s : schemes) {
                            if (s.v > 0){
                                bestSchemes.add(s);
                            }
                        }
                        if (it.overTime) {
                            return true;
                        }
                    }

                    System.out.println("Fract | Lower bound: " + LB
                            + " | Upper bound: " + UB + " | Gap: "
                            + ((UB - LB) / UB) + " | BB Depth: "
                            + depth + " | Local CG: " + CGobj + " | " + "schemes: " + schemes.size());
                    System.out.flush();

                    if (branchValue == 0) {
                        BBtree newNode0 = new BBtree();
                        newNode0.father = bbt;
                        newNode0.son0 = null;
                        newNode0.son1 = null;
                        newNode0.bArcS = selectArc1;
                        newNode0.bArcT = selectArc2;
                        newNode0.branch = branchValue;
                        newNode0.nodeUB = 1E10;
                        newNode0.rBranch = isRB;
                        newNode0.IntV = new int[it.numN];
                        System.arraycopy(bbt.IntV,0,newNode0.IntV,0,it.numN);
                        arcBranching(it, newNode0, true);
                        ArrayList<scheme> node0Schemes = new ArrayList<scheme>();
                        for (scheme s : schemes) {
                            boolean accept = true;
                            if (isRB) {
                                ArrayList<Integer> path = s.r.getpath();
                                iNode = 0;
                                for (i = 1; accept && (i < path.size()); i++) {
                                    jNode = path.get(i);
                                    if ((iNode == selectArc1) && (jNode == selectArc2)) {
                                        accept =false;
                                    }
                                    iNode = jNode;
                                }
                            }
                            else {
                                for (move m : s.mList) {
                                    ArrayList<Integer> group = m.getgroup();
                                    iNode = group.get(0);
                                    for (i = 1; accept && (i < group.size()); i++) {
                                        jNode = group.get(i);
                                        if ((iNode == selectArc2) && (jNode == selectArc1)) {
                                            accept =false;
                                        }
                                    }
                                }
                            }
                            if (accept) {
                                node0Schemes.add(s);
                            }
                        }
                        boolean subTree = BBNode(it,node0Schemes,bestSchemes,newNode0,depth + 1);
                        node0Schemes = null;
                        if (!subTree) {
                            return false;
                        }
                        bbt.son0 = newNode0;

                        BBtree newNode1 = new BBtree();
                        newNode1.father = bbt;
                        newNode1.son0 = null;
                        newNode1.son1 = null;
                        newNode1.bArcS = selectArc1;
                        newNode1.bArcT = selectArc2;
                        newNode1.branch = 1 - branchValue;
                        newNode1.nodeUB =1E10;
                        newNode1.rBranch = isRB;
                        newNode1.IntV = new int[it.numN];
                        System.arraycopy(bbt.IntV,0,newNode1.IntV,0,it.numN);
                        for (i = 0; i < it.numN; i++) {
                            System.arraycopy(it.T[i],0,it.Trb[i],0,it.numN);
                            System.arraycopy(it.T[i],0,it.Tmb[i],0,it.numN);
                        }
                        arcBranching(it, newNode1, true);
                        ArrayList<scheme> node1Schemes = new ArrayList<scheme>();
                        for (scheme s : schemes) {
                            boolean accept = true;
                            if (isRB) {
                                ArrayList<Integer> path = s.r.getpath();
                                iNode = 0;
                                for (i = 1; accept && (i < path.size()); i++) {
                                    jNode = path.get(i);
                                    if (iNode==selectArc1 && jNode != selectArc2) {
                                        accept =false;
                                    } else if (iNode != selectArc1 && jNode == selectArc2) {
                                        accept =false;
                                    }
                                    iNode = jNode;
                                }
                                if (accept) {
                                    for (move m : s.mList) {
                                        ArrayList<Integer> group = m.getgroup();
                                        iNode = group.get(0);
                                        for (i = 1; accept && (i < group.size()); i++) {
                                            jNode = group.get(i);
                                            if (jNode == selectArc1 || jNode == selectArc2) {
                                                accept =false;
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                if (s.r.path.contains(selectArc1)) {
                                    accept =false;
                                }
                                if (accept) {
                                    for (move m : s.mList) {
                                        ArrayList<Integer> group = m.getgroup();
                                        iNode = group.get(0);
                                        for (i = 1; accept && (i < group.size()); i++) {
                                            jNode = group.get(i);
                                            if (jNode == selectArc1 && iNode != selectArc2) {
                                                accept =false;
                                            }
                                        }
                                    }
                                }
                            }
                            if (accept) {
                                node1Schemes.add(s);
                            }
                        }

                        subTree = BBNode(it,node1Schemes,bestSchemes,newNode1,depth + 1);
                        node1Schemes = null;
                        bbt.son1 = newNode1;
                        bbt.nodeUB = Math.max(newNode0.nodeUB, newNode1.nodeUB);
                        return subTree;
                    }

                    else {
                        BBtree newNode1 = new BBtree();
                        newNode1.father = bbt;
                        newNode1.son0 = null;
                        newNode1.son1 = null;
                        newNode1.bArcS = selectArc1;
                        newNode1.bArcT = selectArc2;
                        newNode1.branch = branchValue;
                        newNode1.nodeUB =1E10;
                        newNode1.rBranch = isRB;
                        newNode1.IntV = new int[it.numN];
                        System.arraycopy(bbt.IntV,0,newNode1.IntV,0,it.numN);
                        arcBranching(it, newNode1, true);
                        ArrayList<scheme> node1Schemes = new ArrayList<scheme>();
                        for (scheme s : schemes) {
                            boolean accept = true;
                            if (isRB) {
                                ArrayList<Integer> path = s.r.getpath();
                                iNode = 0;
                                for (i = 1; accept && (i < path.size()); i++) {
                                    jNode = path.get(i);
                                    if (iNode==selectArc1 && jNode != selectArc2) {
                                        accept =false;
                                    } else if (iNode != selectArc1 && jNode == selectArc2) {
                                        accept =false;
                                    }
                                    iNode = jNode;
                                }
                                if (accept) {
                                    for (move m : s.mList) {
                                        ArrayList<Integer> group = m.getgroup();
                                        iNode = group.get(0);
                                        for (i = 1; accept && (i < group.size()); i++) {
                                            jNode = group.get(i);
                                            if (jNode == selectArc1 || jNode == selectArc2) {
                                                accept =false;
                                            }
                                        }
                                    }
                                }
                            }

                            else {
                                if (s.r.path.contains(selectArc1)) {
                                    accept =false;
                                }
                                if (accept) {
                                    for (move m : s.mList) {
                                        ArrayList<Integer> group = m.getgroup();
                                        iNode = group.get(0);
                                        for (i = 1; accept && (i < group.size()); i++) {
                                            jNode = group.get(i);
                                            if (jNode == selectArc1 && iNode != selectArc2) {
                                                accept =false;
                                            }
                                        }
                                    }
                                }
                            }
                            if (accept) {
                                node1Schemes.add(s);
                            }
                        }

                        boolean subTree = BBNode(it,node1Schemes,bestSchemes,newNode1,depth + 1);
                        node1Schemes = null;
                        if (!subTree) {
                            return false;
                        }
                        bbt.son1 = newNode1;

                        BBtree newNode0 = new BBtree();
                        newNode0.father = bbt;
                        newNode0.son0 = null;
                        newNode0.son1 = null;
                        newNode0.bArcS = selectArc1;
                        newNode0.bArcT = selectArc2;
                        newNode0.branch = 1-branchValue;
                        newNode0.nodeUB = 1E10;
                        newNode0.rBranch = isRB;
                        newNode0.IntV = new int[it.numN];
                        System.arraycopy(bbt.IntV,0,newNode0.IntV,0,it.numN);
                        for (i = 0; i < it.numN; i++) {
                            System.arraycopy(it.T[i],0,it.Trb[i],0,it.numN);
                            System.arraycopy(it.T[i],0,it.Tmb[i],0,it.numN);
                        }
                        arcBranching(it, newNode0, true);
                        ArrayList<scheme> node0Schemes = new ArrayList<scheme>();
                        for (scheme s : schemes) {
                            boolean accept = true;
                            if (isRB) {
                                ArrayList<Integer> path = s.r.getpath();
                                iNode = 0;
                                for (i = 1; accept && (i < path.size()); i++) {
                                    jNode = path.get(i);
                                    if ((iNode == selectArc1) && (jNode == selectArc2)) {
                                        accept =false;
                                    }
                                    iNode = jNode;
                                }
                            }
                            else {
                                for (move m : s.mList) {
                                    ArrayList<Integer> group = m.getgroup();
                                    iNode = group.get(0);
                                    for (i = 1; accept && (i < group.size()); i++) {
                                        jNode = group.get(i);
                                        if ((iNode == selectArc2) && (jNode == selectArc1)) {
                                            accept =false;
                                        }
                                    }
                                }
                            }
                            if (accept) {
                                node0Schemes.add(s);
                            }
                        }

                        subTree = BBNode(it,node0Schemes,bestSchemes,newNode0,depth + 1);
                        node0Schemes = null;
                        bbt.son0 = newNode0;
                        bbt.nodeUB = Math.max(newNode0.nodeUB, newNode1.nodeUB);
                        return subTree;
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error: " + e);
        }
        return false;
    }
}
