package AlgorithmicGroupTheory;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class connectivityTest {
    /** 3.6.2. Connectivity Test */

    /**
     * Finding the N values as explained in 3.6.2.
     *
     * @param index int row index
     * @param total int number of atoms.
     * @param mat int[][] adjacency matrix
     * @return Set<Integer>
     */
    public static Set<Integer> nValues(int index, int total, int[][] mat) {
        Set<Integer> nValues = new HashSet<Integer>();
        nValues.add(index);
        int[] theRow = mat[index];
        for (int i = (index + 1); i < total; i++) {
            if (theRow[i] > 0) {
                nValues.add(i);
            }
        }
        return nValues;
    }

    /**
     * Finding the W values as explained in 3.6.2.
     *
     * @param nValues ArrayList<Integer> N values
     * @param Kformer ArrayList<Integer> the K values of the former step
     * @return Set<Integer>
     */
    public static Set<Integer> wValues(Set<Integer> nValues, int[] Kformer) {
        Set<Integer> wValues = new HashSet<Integer>();
        for (Integer i : nValues) {
            wValues.add(Kformer[i]);
        }
        return wValues;
    }

    /**
     * Finding the K values as explained in 3.6.2.
     *
     * @param wValues Set<Integer> wValues
     * @param kFormer ArrayList<Integer> the K values of the former step
     * @return ArrayList<Integer>
     */
    public static int[] kValues(int total, Set<Integer> wValues, int[] kFormer) {
        int[] kValues = new int[total];
        int min = Collections.min(wValues);
        for (int i = 0; i < total; i++) {
            if (wValues.contains(kFormer[i])) {
                kValues[i] = min;
            } else {
                kValues[i] = kFormer[i];
            }
        }
        return kValues;
    }

    /**
     * Initializing the first K values list.
     *
     * @param total int number of atoms.
     * @return ArrayList<Integer>
     */
    public static int[] initialKList(int total) {
        int[] k = new int[total];
        for (int i = 0; i < total; i++) {
            k[i] = i;
        }
        return k;
    }

    /**
     * Test whether an adjacency matrix is connected or disconnected.
     *
     * @param p int the index where the atoms with degree 1 starts.
     * @param mat int[][] adjacency matrix
     * @return boolean
     */
    public static boolean connectivityTest(int p, int[][] mat) {
        boolean check = false;
        int total = mat.length;
        int[] kValues = initialKList(total);
        Set<Integer> nValues = new HashSet<Integer>();
        Set<Integer> wValues = new HashSet<Integer>();
        int zValue = 0;
        for (int i = 0; i < p; i++) {
            nValues = nValues(i, total, mat);
            wValues = wValues(nValues, kValues);
            zValue = Collections.min(wValues);
            kValues = kValues(total, wValues, kValues);
        }
        if (zValue == 0 && allIs0(kValues)) {
            check = true;
        }
        return check;
    }

    /**
     * Checks whether all the entries are equal to 0 or not. (Grund 3.6.2)
     *
     * @param list ArrayList<Integer>
     * @return boolean
     */
    public static boolean allIs0(int[] list) {
        boolean check = true;
        for (int i = 0; i < list.length; i++) {
            if (list[i] != 0) {
                check = false;
                break;
            }
        }
        return check;
    }
}
