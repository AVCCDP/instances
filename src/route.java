import java.util.ArrayList;

public class route implements Cloneable {
    public double profitR;
    public double inspectR;
    public double rQ;
    public double v;
    public ArrayList<Integer> path;
    public ArrayList<Double> pathAT;
    public ArrayList<Double> pathLT;

    
    public route() {
        this.path = new ArrayList<Integer>();
        this.profitR = 0.0;
        this.pathAT = new ArrayList<Double>();
        this.pathLT = new ArrayList<Double>();
    }
    
    public double getProfit(Instance it) {
        this.profitR=0.0;
        int i,jNode;

        for (i = 1; i < this.path.size()-1; i++) {
            jNode = this.path.get(i);
            this.profitR += it.P[jNode]-this.pathAT.get(i)*it.d[jNode];
        }
        return this.profitR;
    }
    
    public double getInspectR(Instance it) {
        int i,jNode;
        this.inspectR = it.ps;

        for (i = 1; i < this.path.size()-1; i++) {
            jNode = this.path.get(i);
            this.inspectR += it.P_r[jNode] + this.pathAT.get(i)*it.d[jNode];
        }
        return this.inspectR;
    }
    
    public void updateLT(int ind, double LT, Instance it,scheme newS) {
        this.pathLT.set(ind,LT);
        ind++;
        for (int i = ind; i < this.path.size(); i++) {
            int prei = this.path.get(i-1);
            int curi = this.path.get(i);
            double newAT = this.pathLT.get(i-1)+it.T[prei][curi];
            this.pathAT.set(i,newAT);

            if (this.pathAT.get(i) > this.pathLT.get(i)) {
                this.pathLT.set(i,newAT);
            }

        }

        for (move m : newS.mList) {
            int meet = m.getgroup().get(0);
            int meetInd = this.path.indexOf(meet);

            if (meetInd >= ind) {
                for (int i = 1; i < m.group.size(); i++) {

                    if (m.groupAT.get(i) < this.pathAT.get(meetInd)) {
                        m.getgroupAT().set(i,this.pathAT.get(meetInd));
                    }
                }
            }
        }
    }
    
    public ArrayList<Integer> getpath() {
        return this.path;
    }
    
    public void addcity(int city, Instance it) {
        this.path.add(city);

        if (city == 0) {
            this.pathAT.add(0.0);
            this.pathLT.add(0.0);
        } else {
            int prei = this.path.get(this.path.size()-2);
            double Time = Math.round((this.pathLT.get(this.path.size()-2)+it.T[prei][city])*100)/100.0;
            this.pathAT.add(Time);
            this.pathLT.add(Time);
        }
    }
    
    public void setInspectR(double c) {
        this.inspectR = c;
    }
    
    public void setV(double a) {
        this.v = a;
    }
    
    public double getV() {
        return this.v;
    }
    
    public ArrayList<Double> getpathAT() {
        return this.pathAT;
    }
    
    public ArrayList<Double> getpathLT() {
        return this.pathLT;
    }
    
    public void setRQ(double q) {
        this.rQ = q;
    }
    
    public route clone() throws CloneNotSupportedException {
        route route = (route) super.clone();
        route.profitR = this.profitR;
        route.inspectR = this.inspectR;
        route.path = (ArrayList<Integer>) this.path.clone();
        route.pathAT = (ArrayList<Double>) this.pathAT.clone();
        route.pathLT = (ArrayList<Double>) this.pathLT.clone();
        return route;
    }
}
