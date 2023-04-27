package org.xcgis.hpgc.isodata;

import org.xcgis.hpgc.Tuple;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public class ParallelISODATA {
    private double[][] data;
    private List<double[]> centers;
    private List<List<Integer>> clusters;
    private int K, thetaN, thetaS, thetaC, L, I, dim, partition, partitionNum;
    private Lock lock = new ReentrantLock();


    public ParallelISODATA(double[][] data, List<double[]> centers,
                   int K, int thetaN, int thetaS, int thetaC, int L, int I, int partition) {
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
        this.partition = partition;
        this.partitionNum = data.length / partition;
    }

    public int[] getResult() throws InterruptedException {
        process();

        int[] result = new int[data.length];
        for (int j = 0; j < centers.size(); ++j) {
            for (int index : clusters.get(j)) {
                result[index] = j + 1;
            }
        }

        return result;
    }

    public void process() throws InterruptedException {
        for (int count = 0; count < I; ++count) {
            System.out.println("第" + (count + 1) + "次聚类");
            clump();
            checkClusterCount();
            updateCenters();
            System.out.print("第" + (count + 1) + "轮循环为");
            if (centers.size() <= K / 2 || (count % 2 == 0 && centers.size() < 2 * K)) {
                System.out.println("分裂操作");
                divide();
            } else if (clusters.size() > 2 * K || count % 2 == 1) {
                System.out.println("合并操作");
                merge();
            }
        }
    }

    public void clump() throws InterruptedException {
        for (List<Integer> cluster : clusters) {
            cluster.clear();
        }

        CountDownLatch countDownLatch = new CountDownLatch(partition);
        for (int p = 0; p < partition; ++p) {
            int finalP = p;
            new Thread(() -> {
                int end = finalP == partition - 1 ? data.length : (finalP + 1) * partitionNum;
                for (int i = finalP * partitionNum; i < end; ++i) {
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
                    lock.lock();
                    clusters.get(clusterIndex).add(i);
                    lock.unlock();
                }
                countDownLatch.countDown();
            }).start();
        }
        countDownLatch.await();
    }

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

    public void updateCenters() throws InterruptedException {
        CountDownLatch countDownLatch = new CountDownLatch(centers.size());
        for (int j = 0; j < centers.size(); ++j) {
            int finalJ = j;
            new Thread(() -> {
                List<Integer> cluster = clusters.get(finalJ);
                for (int k = 0; k < dim; ++k) {
                    double sum = 0;
                    for (int i = 0; i < cluster.size(); ++i) {
                        sum += data[cluster.get(i)][k];
                    }
                    centers.get(finalJ)[k] = sum / cluster.size();
                }
                countDownLatch.countDown();
            }).start();
        }
        countDownLatch.await();

//        int m = 1;
//        for (List<Integer> cluster : clusters) {
//            System.out.print("第" + m + "类(");
//            for (int index : cluster) {
//                System.out.print("x" + (index + 1) + ",");
//            }
//            System.out.println(")");
//            m += 1;
//        }
    }

    public Tuple<double[], Double> getAverageDistance() throws InterruptedException {
        AtomicReference<Double> averageDistance = new AtomicReference<>((double) 0);
        double[] adInClusters = new double[centers.size()];

        CountDownLatch countDownLatch = new CountDownLatch(centers.size());
        for (int j = 0; j < centers.size(); ++j) {
            List<Integer> cluster = clusters.get(j);
            int finalJ = j;
            new Thread(() -> {
                double adInCluster = 0;
                for (int i = 0; i < cluster.size(); ++i) {
                    double adSingle = 0;
                    for (int k = 0; k < dim; ++k) {
                        adSingle += Math.pow(data[cluster.get(i)][k] - centers.get(finalJ)[k], 2);
                    }
                    adInCluster += Math.pow(adSingle, 0.5);
                }
                adInClusters[finalJ] = adInCluster;
                double finalAdInCluster = adInCluster;
                averageDistance.updateAndGet(v -> new Double((double) (v + finalAdInCluster))); // 原子操作
                countDownLatch.countDown();
            }).start();
        }
        countDownLatch.await();

        return new Tuple<>(adInClusters, averageDistance.get() / data.length);
    }

    public Tuple<double[], int[]> getMaxSDVectorComponents() throws InterruptedException {
        double[][] sdVectors = new double[centers.size()][dim];
        double[] maxSDVectorComponents = new double[centers.size()];
        int[] maxSDVectorComponentsIndex = new int[centers.size()];

        CountDownLatch countDownLatch = new CountDownLatch(centers.size());
        for (int j = 0; j < centers.size(); ++j) {
            List<Integer> cluster = clusters.get(j);
            double[] center = centers.get(j);
            int finalJ = j;
            new Thread(() -> {
                for (int k = 0; k < dim; ++k) {
                    for (int i = 0; i < cluster.size(); ++i) {
                        sdVectors[finalJ][k] += Math.pow(data[cluster.get(i)][k] - center[k], 2);
                    }
                    sdVectors[finalJ][k] = Math.pow(sdVectors[finalJ][k] / cluster.size(), 0.5);
                    if (sdVectors[finalJ][k] > maxSDVectorComponents[finalJ]) {
                        maxSDVectorComponents[finalJ] = sdVectors[finalJ][k];
                        maxSDVectorComponentsIndex[finalJ] = k;
                    }
                }
                countDownLatch.countDown();
            }).start();
        }
        countDownLatch.await();

        return new Tuple<>(maxSDVectorComponents, maxSDVectorComponentsIndex);
    }

    public void divide() throws InterruptedException {
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
    }

    public void merge() {
        System.out.println("合并前的中心点如下：");
        for (double[] center : centers) {
            System.out.print("(");
            for (double v : center) {
                System.out.print(v + ",");
            }
            System.out.println(")");
        }

        int oldSize = centers.size();
        for (int i = 0; i < oldSize; ++i) {
            for (int j = i + 1; j < centers.size(); ++j) {
                double distance = 0;
                for (int k = 0; k < dim; ++k) {
                    distance += Math.pow(centers.get(i)[k] - centers.get(j)[k], 2);
                }
                distance = Math.pow(distance, 0.5);
                if (distance < thetaC) {
                    int ni = clusters.get(i).size();
                    int nj = clusters.get(j).size();
                    for (int k = 0; k < dim; ++k) {
                        centers.get(i)[k] = (ni * centers.get(i)[k] + nj * centers.get(j)[k]) / (ni + nj);
                    }
                    centers.remove(j);
                    clusters.get(i).addAll(clusters.get(j));
                    clusters.remove(j);
                    --j;
                }
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
