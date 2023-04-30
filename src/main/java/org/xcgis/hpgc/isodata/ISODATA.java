package org.xcgis.hpgc.isodata;

import org.xcgis.hpgc.Main;
import org.xcgis.hpgc.Tuple;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * 串行版ISODATA算法，在本类的方法中，循环里的j一般为centers的index，k一般为波段数的循环index
 */
public class ISODATA {
    // 影像数据
    private double[][] data;
    // 聚类中心
    private List<double[]> centers;
    // 聚类结果，只记录被聚类点的索引
    private List<List<Integer>> clusters;
    // 分别为：预期的聚类中心数目，每一聚类域中最少的样本数目，一个聚类域中样本距离分布的标准差
    // 两个聚类中心间的最小距离，在一次迭代运算中可以合并的聚类中心的最多对数，迭代运算的次数，波段数
    private int K, thetaN, thetaS, thetaC, L, I, dim;

    public ISODATA(double[][] data, List<double[]> centers,
                   int K, int thetaN, int thetaS, int thetaC, int L, int I) {
        this.data = data;
        this.centers = centers;
        this.clusters = new ArrayList<>();
        for (int i = 0; i < centers.size(); ++i) {
            clusters.add(new ArrayList<>());
        }
        this.K = K;
        this.thetaN = thetaN;
        this.thetaS = thetaS;
        this.thetaC = thetaC;
        this.L = L;
        this.I = I;
        this.dim = data[0].length;
    }

    /**
     * 获取聚类结果，按照clusters中记录的索引组织result
     *
     * @return 最终聚类结果
     */
    public int[] getResult() {
        process();

        int[] result = new int[data.length];
        for (int j = 0; j < centers.size(); ++j) {
            for (int index : clusters.get(j)) {
                result[index] = j + 1;
            }
        }

        return result;
    }

    /**
     * 迭代计算过程
     */
    public void process() {
        for (int count = 0; count < I; ++count) {
            System.out.println("第" + (count + 1) + "次聚类");
            clump();
            checkClusterCount();
            updateCenters();
            System.out.print("第" + (count + 1) + "轮循环为");
            if (centers.size() <= K / 2 || (count % 2 == 0 && centers.size() < 2 * K)) {
                System.out.println("分裂操作");
                int check = divide();
                if (check == 0) {
                    System.out.println("分裂未完成，开始合并");
                    merge();
                }
            } else if (clusters.size() > 2 * K || count % 2 == 1) {
                System.out.println("合并操作");
                merge();
            }
        }
    }

    /**
     * 将所有模式样本分给最近的聚类
     */
    public void clump() {
        for (List<Integer> cluster : clusters) {
            cluster.clear();
        }

        for (int i = 0; i < data.length; ++i) {
            double minDis = Double.MAX_VALUE;
            int clusterIndex = 0;
            for (int j = 0; j < centers.size(); ++j) {
                double sum = 0;
                for (int k = 0; k < dim; ++k) {
                    sum += Math.pow(data[i][k] - centers.get(j)[k], 2);
                }
                sum = Math.pow(sum, 0.5);
                if (sum < minDis) {
                    minDis = sum;
                    clusterIndex = j;
                }
            }
            clusters.get(clusterIndex).add(i);
        }
    }

    /**
     * 如果clusters(j)中的样本数目小于thetaN，则取消该样本子集
     */
    public void checkClusterCount() {
        for (int i = 0; i < centers.size(); ++i) {
            if (clusters.get(i).size() < thetaN) {
                for (int index : clusters.get(i)) {
                    double minDis = Double.MAX_VALUE;
                    int clusterIndex = 0;
                    for (int j = 0; j < centers.size(); ++j) {
                        if (j != i) {
                            double sum = 0;
                            for (int k = 0; k < dim; ++k) {
                                sum += Math.pow(data[index][k] - centers.get(j)[k], 2);
                            }
                            sum = Math.pow(sum, 0.5);
                            if (sum < minDis) {
                                minDis = sum;
                                clusterIndex = j;
                            }
                        }
                    }
                    clusters.get(clusterIndex).add(index);
                }
                centers.remove(i);
                clusters.remove(i);
                --i;
            }
        }
    }

    /**
     * 修正聚类中心
     */
    public void updateCenters() {
        for (int j = 0; j < centers.size(); ++j) {
            List<Integer> cluster = clusters.get(j);
            for (int k = 0; k < dim; ++k) {
                double sum = 0;
                for (int i = 0; i < cluster.size(); ++i) {
                    sum += data[cluster.get(i)][k];
                }
                centers.get(j)[k] = sum / cluster.size();
            }
        }
    }

    /**
     * 计算各聚类域中模式样本与各聚类中心间的平均距离（此步骤中用平均距离与聚类数目的乘积做代替，即距离总和）
     * 及全部模式样本和其对应聚类中心的总平均距离
     *
     * @return 各聚类域中模式样本与各聚类中心间的距离总和与全部模式样本和其对应聚类中心的总平均距离的二元组
     */
    public Tuple<double[], Double> getAverageDistance() {
        double averageDistance = 0;
        double[] adInClusters = new double[centers.size()];

        for (int j = 0; j < centers.size(); ++j) {
            double adInCluster = 0;
            for (int i = 0; i < clusters.get(j).size(); ++i) {
                double adSingle = 0;
                for (int k = 0; k < dim; ++k) {
                    adSingle += Math.pow(data[clusters.get(j).get(i)][k] - centers.get(j)[k], 2);
                }
                adInCluster += Math.pow(adSingle, 0.5);
            }
            adInClusters[j] = adInCluster;
            averageDistance += adInCluster;
        }
        averageDistance /= data.length;
        return new Tuple<>(adInClusters, averageDistance);
    }

    /**
     * 计算每个聚类中样本距离的标准差向量，且求出每一标准差向量中的最大分量
     *
     * @return 各聚类标准差向量的最大分量及它们对应的波段数的数组二元组
     */
    public Tuple<double[], int[]> getMaxSDVectorComponents() {
        double[][] sdVectors = new double[centers.size()][dim];
        double[] maxSDVectorComponents = new double[centers.size()];
        int[] maxSDVectorComponentsIndex = new int[centers.size()];

        for (int j = 0; j < centers.size(); ++j) {
            for (int k = 0; k < dim; ++k) {
                for (int i = 0; i < clusters.get(j).size(); ++i) {
                    sdVectors[j][k] += Math.pow(data[clusters.get(j).get(i)][k] - centers.get(j)[k], 2);
                }
                sdVectors[j][k] = Math.pow(sdVectors[j][k] / clusters.get(j).size(), 0.5);
                if (sdVectors[j][k] > maxSDVectorComponents[j]) {
                    maxSDVectorComponents[j] = sdVectors[j][k];
                    maxSDVectorComponentsIndex[j] = k;
                }
            }
        }

        return new Tuple<>(maxSDVectorComponents, maxSDVectorComponentsIndex);
    }

    /**
     * 聚类的分裂操作
     *
     * @return 用于判断分裂操作是否完成的标记，若结果为0，说明分裂操作未完成
     */
    public int divide() {
        System.out.println("分裂前的中心点如下：");
        for (double[] center : centers) {
            System.out.print("(");
            for (double v : center) {
                System.out.print(v + ",");
            }
            System.out.println(")");
        }

        Tuple<double[], Double> t1 = getAverageDistance();
        double[] adInClusters = t1.getV1();
        double averageDistance = t1.getV2();
        Tuple<double[], int[]> t2 = getMaxSDVectorComponents();
        double[] maxSDVectorComponents = t2.getV1();
        int[] maxSDVectorComponentsIndex = t2.getV2();

        int check = 0;
        int oldSize = centers.size();
        for (int j = 0; j < oldSize; ++j) {
            List<Integer> cluster = clusters.get(j);
            if (maxSDVectorComponents[j] > thetaS && (
                (cluster.size() > 2 * thetaN && adInClusters[j] / cluster.size() > averageDistance) ||
                centers.size() <= K / 2
            )) {
                double[] newCenter = new double[dim];
                for (int k = 0; k < dim; ++k) {
                    if (k == maxSDVectorComponentsIndex[j]) {
                        newCenter[k] = centers.get(j)[k] - 0.5 * maxSDVectorComponents[j];
                        centers.get(j)[k] += 0.5 * maxSDVectorComponents[j];
                    } else {
                        newCenter[k] = centers.get(j)[k];
                    }
                }
                centers.add(newCenter);
                clusters.add(new ArrayList<>());
                ++check;
            }
        }

        System.out.println("分裂后的中心点如下：");
        for (double[] center : centers) {
            System.out.print("(");
            for (double v : center) {
                System.out.print(v + ",");
            }
            System.out.println(")");
        }

        return check;
    }

    /**
     * 聚类的合并操作
     */
    public void merge() {
        System.out.println("合并前的中心点如下：");
        for (double[] center : centers) {
            System.out.print("(");
            for (double v : center) {
                System.out.print(v + ",");
            }
            System.out.println(")");
        }

        List<Tuple<Tuple<Integer, Integer>, Double>> distancesWithIndex = new ArrayList<>();
        for (int i = 0; i < centers.size(); ++i) {
            for (int j = i + 1; j < centers.size(); ++j) {
                double distance = 0;
                for (int k = 0; k < dim; ++k) {
                    distance += Math.pow(centers.get(i)[k] - centers.get(j)[k], 2);
                }
                distance = Math.pow(distance, 0.5);
                if (distance < thetaC) {
                    distancesWithIndex.add(new Tuple<>(new Tuple<>(i, j), distance));
                }
            }
        }

        distancesWithIndex.sort(Comparator.comparing(Tuple::getV2));
        int l = 0;
        int[] check = new int[centers.size()];
        for (int m = 0; m < distancesWithIndex.size(); ++m) {
            if (l == L) {
                break;
            }
            int i = distancesWithIndex.get(m).getV1().getV1();
            int j = distancesWithIndex.get(m).getV1().getV2();
            if (check[i] != 1 && check[j] != 1) {
                int ni = clusters.get(i).size();
                int nj = clusters.get(j).size();
                for (int k = 0; k < dim; ++k) {
                    centers.get(i)[k] = (ni * centers.get(i)[k] + nj * centers.get(j)[k]) / (ni + nj);
                }
                centers.remove(j);
                clusters.get(i).addAll(clusters.get(j));
                clusters.remove(j);
                check[i] = 1;
                check[j] = 1;
                ++l;
            }
        }

        System.out.println("合并后的中心点如下：");
        for (double[] center : centers) {
            System.out.print("(");
            for (double v : center) {
                System.out.print(v + ",");
            }
            System.out.println(")");
        }
    }
}
