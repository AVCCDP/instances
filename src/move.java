import java.util.ArrayList;

public class move implements Cloneable {
    public double profitM;
    public double inspectM;
    public double v;
    public ArrayList<Integer> group;
    public ArrayList<Double> groupAT;
    
    public move() {
        this.profitM = 0.0;
        this.group = new ArrayList<Integer>();
        this.groupAT = new ArrayList<Double>();
    }
    
    public void addcity (int city,Instance it) {
        this.group.add(city);
        if (this.group.size() == 1) {
            this.groupAT.add(0.0);
        } else {
            this.groupAT.add(it.T[city][this.group.get(0)]);
        }
    }
    
    public double getProfit(Instance it) {
        this.profitM = 0.0;
        int i, mNode, jNode;
        mNode = this.group.get(0);
        for (i = 1; i < this.group.size(); i++) {
            jNode = this.group.get(i);
            this.profitM += it.P[jNode]*it.md[jNode][mNode]-this.groupAT.get(i)*it.d[jNode];
        }
        return this.profitM;
    }
    
    public double getInspectM(Instance it) {
        this.inspectM = 0.0;
        int i, mNode, jNode;
        mNode = this.group.get(0);
        for (i = 1; i < this.group.size(); i++) {
            jNode = this.group.get(i);
            this.inspectM += it.P_m[jNode] - it.P[jNode]*it.md[jNode][mNode] + this.groupAT.get(i)*it.d[jNode];
        }
        return inspectM;
    }
    
    public ArrayList<Integer> getgroup() {
        return this.group;
    }
    
    public void updateAT(route r, Instance it, scheme newS) {
        int mind = r.getpath().indexOf(this.group.get(0));
        for (int i = 1; i < this.group.size(); i++) {
            if (this.groupAT.get(i) > r.pathLT.get(mind)) {
                r.updateLT(mind,this.groupAT.get(i),it,newS);
            }
            else if (this.groupAT.get(i) < r.pathAT.get(mind)) {
                this.groupAT.set(i,r.pathAT.get(mind));
            }
        }
    }
    
    public void setV(double a) {
        this.v = a;
    }
    
    public double getV() {
        return this.v;
    }
    
    public ArrayList<Double> getgroupAT() {
        return this.groupAT;
    }
    
    public move clone() throws CloneNotSupportedException {
        move move = (move) super.clone();
        move.profitM = this.profitM;
        move.inspectM = this.inspectM;
        move.group = (ArrayList<Integer>) this.group.clone();
        move.groupAT = (ArrayList<Double>) this.groupAT.clone();
        return move;
    }

}
