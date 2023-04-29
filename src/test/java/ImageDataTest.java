import org.xcgis.hpgc.io.ImageIO;
import org.xcgis.hpgc.isodata.ISODATA;

import java.util.ArrayList;
import java.util.List;

public class ImageDataTest {
    public static void main(String[] args) {
        ImageIO imageIO = new ImageIO();
        double[][] data = imageIO.read("D:\\Projects\\Java\\python-java\\band_clip_6.tif");
        List<double[]> centers = new ArrayList<>();
        centers.add(new double[]{0, 0, 0, 0, 0});
        ISODATA isodata = new ISODATA(data, centers, 6, data.length / 24, 1, 5, 1, 10);
        long start = System.currentTimeMillis();
        int[] ans = isodata.getResult();
        long end = System.currentTimeMillis();
        System.out.println("用时：" + (end - start));
        // imageIO.write(ans, "D:\\xx1.tif");
    }
}
