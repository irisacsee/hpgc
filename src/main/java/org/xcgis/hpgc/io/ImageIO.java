package org.xcgis.hpgc.io;

import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.Driver;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;

/**
 * 影像读取写入操作封装类
 */
public class ImageIO {
    private int width;
    private int height;
    private String projection;
    private double[] transform;

    public ImageIO() {}

    public double[][] read(String filePath) {
        gdal.AllRegister();
        Dataset source = gdal.Open(filePath, gdalconst.GA_ReadOnly);

        int bands = source.getRasterCount();
        width = source.getRasterXSize();
        height = source.getRasterYSize();
        projection = source.GetProjection();
        transform = source.GetGeoTransform();
        System.out.println("width: " + width + ", height: " + height);

        double[][] bandValues = new double[bands][width * height];
        for (int k = 0; k < bands; ++k) {
            source.GetRasterBand(k + 1).ReadRaster(0, 0, width, height, bandValues[k]);
        }
        source.delete();

        double[][] ans = new double[width * height][bands];
        for (int i = 0; i < width * height; ++i) {
            for (int k = 0; k < bands; ++k) {
                ans[i][k] = bandValues[k][i];
            }
        }

        return ans;
    }

    public void write(int[] result, String filePath) {
        Dataset target = gdal.GetDriverByName("GTiff").Create(filePath, width, height, 1);
        target.SetProjection(projection);
        target.SetGeoTransform(transform);
        Band band = target.GetRasterBand(1);
        band.WriteRaster(0, 0, width, height, result);
        band.FlushCache();
        target.delete();
    }
}
