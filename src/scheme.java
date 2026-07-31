import java.util.ArrayList;

public class scheme {
    public double profitS;
    public double inspectS;
    public double Q;
    public double v;
    public ArrayList<Integer> customers;
    public route r;
    public ArrayList<move> mList;

    
    public scheme() {
        this.r = new route();
        this.mList = new ArrayList<move>();
        this.customers = new ArrayList<Integer>();
    }
    
    public void setR(route r) {
        this.r = r;
        for (int i = 1; i < r.path.size()-1; i++) {
            this.customers.add(r.path.get(i));
        }
    }
    
    public void setM(move m) {
        this.mList.add(m);
        for (int i = 1; i < m.group.size(); i++) {
            this.customers.add(m.group.get(i));
        }
    }
    
    public void setQ(double q){this.Q = q;}
    
    public double getProfit(Instance it) {
        this.profitS =0.0;
        for (move m : mList) {
            this.profitS += m.getProfit(it);
        }
        this.profitS += r.getProfit(it);
        return this.profitS;
    }
    
    public void setInspectS(double inspectS) {
        this.inspectS = inspectS;
    }
    
    public void setV(double a) {
        this.v = a;
    }
    
    public double getV() {
        return this.v;
    }
    
    public void updateCust(){
        for (move m : this.mList) {
            for (int i = 1; i < m.group.size(); i++) {
                if (!this.customers.contains(m.group.get(i))) {
                    this.customers.add(m.group.get(i));
                }
            }
        }
    }
    
    public scheme clone() throws CloneNotSupportedException {
        scheme newScheme = new scheme();
        route newRoute = this.r.clone();
        for (move m : this.mList) {
            newScheme.mList.add(m.clone());
        }
        newScheme.setR(newRoute);
        newScheme.customers = (ArrayList<Integer>) this.customers.clone();
        newScheme.setQ(this.Q);
        return newScheme;
    }

}
