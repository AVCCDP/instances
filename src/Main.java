import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class Main {
    public static void main(String[] args) throws IOException {

        Instance it = new Instance();
        it.readData("input your data file path");
        InitialSolution initSol = new InitialSolution(it);
        ArrayList<scheme> initScheme = initSol.getScheme();
        ArrayList<scheme> bestScheme = new ArrayList<scheme>();
        BranchPrice bp = new BranchPrice();
        bp.BBNode(it,initScheme,bestScheme,null,0);

        double optR = 0.0;
        double optM = 0.0;
        int k = 1;
        for (scheme s : bestScheme) {
            double[] profitR = new double[s.r.path.size()];
            for (int i = 1; i < s.r.path.size()-1; i++) {
                int node = s.r.path.get(i);
                profitR[i] = Math.round((it.P[node] - s.r.pathAT.get(i)*it.d[node])*100)/100.0;
            }
            optR += s.r.profitR;
            if (!s.mList.isEmpty()) {
                for (move m : s.mList) {
                    double[] profitM = new double[m.group.size()];
                    int meet = m.group.get(0);
                    for (int i = 1; i < m.group.size(); i++) {
                        int node = m.group.get(i);
                        profitM[i] = Math.round((it.P[node]*it.md[node][meet] - m.groupAT.get(i)*it.d[node])*100)/100.0;
                    }
                    optM += m.profitM;
                }
            }
        }

        System.out.println("best profit = " + (optR+optM));

    }
}
