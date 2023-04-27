import org.xcgis.hpgc.isodata.ISODATA;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;

public class SmallDataTest {
    public static void main(String[] args) {
        double[][] data = {
                {0, 0}, {3, 8}, {2, 2}, {1, 1}, {5, 3},
                {4, 8}, {6, 3}, {5, 4}, {6, 4}, {7, 5}
        };
        List<double[]> centers = new ArrayList<>();
        centers.add(new double[]{0, 0});
        ISODATA isodata = new ISODATA(data, centers, 3, 1, 1, 4, 1, 4);
        long start = System.currentTimeMillis();
        int[] ans = isodata.getResult();
        long end = System.currentTimeMillis();
        for (int a : ans) {
            System.out.print(a + ", ");
        }
        System.out.println("用时：" + (end - start));
    }
}
