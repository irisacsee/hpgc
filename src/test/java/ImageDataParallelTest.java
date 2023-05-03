import org.xcgis.hpgc.io.ImageIO;
import org.xcgis.hpgc.isodata.ParallelISODATA;

import java.util.ArrayList;
import java.util.List;

public class ImageDataParallelTest {
    public static void main(String[] args) throws InterruptedException {
        ImageIO imageIO = new ImageIO();
        double[][] data = imageIO.read("src/main/resources/band_clip_6.tif");
        List<double[]> centers = new ArrayList<>();
        centers.add(new double[]{0, 0, 0, 0, 0});
        ParallelISODATA parallelISODATA = new ParallelISODATA(
                data, centers, 6, data.length / 24, 1, 5, 1, 10, 5
        );
        long start = System.currentTimeMillis();
        int[] ans = parallelISODATA.getResult();
        long end = System.currentTimeMillis();
        System.out.println("用时：" + (end - start));
        imageIO.write(ans, "src/main/resources/result.tif");
    }
}
