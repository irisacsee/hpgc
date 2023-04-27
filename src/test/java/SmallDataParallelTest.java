import org.xcgis.hpgc.isodata.ISODATA;
import org.xcgis.hpgc.isodata.ParallelISODATA;

import java.util.ArrayList;
import java.util.List;

public class SmallDataParallelTest {
    public static void main(String[] args) throws InterruptedException {
        double[][] data = {
                {0, 0}, {3, 8}, {2, 2}, {1, 1}, {5, 3},
                {4, 8}, {6, 3}, {5, 4}, {6, 4}, {7, 5}
        };
        List<double[]> centers = new ArrayList<>();
        centers.add(new double[]{0, 0});
        ParallelISODATA parallelISODATA = new ParallelISODATA(data, centers, 3, 1, 1, 4, 1, 4, 5);
        long start = System.currentTimeMillis();
        int[] ans = parallelISODATA.getResult();
        long end = System.currentTimeMillis();
        for (int a : ans) {
            System.out.print(a + ", ");
        }
        System.out.println("用时：" + (end - start));
    }
}
