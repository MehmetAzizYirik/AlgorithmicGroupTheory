package AlgorithmicGroupTheory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;

public class maximalTest {
    public List<List<Permutation>> formerPermutations = new ArrayList<List<Permutation>>();
    public List<ArrayList<Integer>> partitionList = new ArrayList<ArrayList<Integer>>();
    public boolean biggest = true;
    public boolean equalAuto = true;

    /** Canonical Test Functions */
    public int findY(int r) {
        return (sum(trialAbstractClass.initialPartition, (r - 1)));
    }

    public int findZ(int r) {
        return (sum(trialAbstractClass.initialPartition, r) - 1);
    }

    /** add number of 1s into a ArrayList */
    public ArrayList<Integer> addOnes(ArrayList<Integer> list, int number) {
        for (int i = 0; i < number; i++) {
            list.add(1);
        }
        return list;
    }

    public boolean canonicalTest(int[][] matrix)
            throws IOException, CloneNotSupportedException, CDKException {
        if (blockTest(trialAbstractClass.r, matrix)) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * To get the canonical partition like in Grund Thesis 3.3.11
     *
     * @param i int row index
     * @param partition ArrayList<Integer> partition
     * @return
     */
    public ArrayList<Integer> canonicalPartition(int i, ArrayList<Integer> partition) {
        return partitionCriteria(partition, i + 1);
    }

    /**
     * Grund Thesis 3.3.2 Partitioning criteria (DONE)
     *
     * @param partEx the former partition
     * @param degree degree of the partitioning.
     * @return
     */
    public ArrayList<Integer> partitionCriteria(ArrayList<Integer> partEx, int degree) {
        ArrayList<Integer> partNew = new ArrayList<Integer>();
        partNew = addOnes(partNew, degree);
        /** I had (degree-1) */
        if (partEx.get(degree - 1) > 1) {
            partNew.add(partEx.get(degree - 1) - 1);
            for (int k = degree; k < partEx.size(); k++) {
                partNew.add(partEx.get(k));
            }
        } else if (partEx.get(degree - 1) == 1) {
            for (int k = degree; k < partEx.size(); k++) {
                partNew.add(partEx.get(k));
            }
        }
        return partNew;
    }

    /**
     * Grund Thesis 3.3.3. To calculate the number of conjugacy classes, used in cycle transposition
     * calculation.
     *
     * @param partEx ArrayList<Integer> former atom partition
     * @param degree ArrayList<Integer> atom valences
     * @return
     */
    public int LValue(ArrayList<Integer> partEx, int degree) {
        return (sum(partEx, (degree))
                - (degree)); // In Grund, the numeration starts with 1. Since we start with 0, it
        // should be -(degree+1)+1, so just -degree
    }

    public List<Permutation> cycleTranspositions(int index, ArrayList<Integer> partition) {
        int total = sum(partition);
        List<Permutation> perms = new ArrayList<Permutation>();
        int lValue = LValue(partition, index);
        for (int i = 0; i < lValue; i++) {
            int[] values = idValues(total);
            int former = values[index];
            values[index] = values[index + i];
            values[index + i] = former;
            Permutation p = new Permutation(values);
            perms.add(p);
        }
        return perms;
    }

    /**
     * Values for an id permutation for a given size
     *
     * @param size
     * @return
     */
    public int[] idValues(int size) {
        int[] id = new int[size];
        for (int i = 0; i < size; i++) {
            id[i] = i;
        }
        return id;
    }

    public boolean allis1(ArrayList<Integer> partition) {
        boolean check = true;
        for (int i = 0; i < partition.size(); i++) {
            if (partition.get(i) != 1) {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * Updating canonical partition list.
     *
     * @param index row index
     * @param newPartition atom partition
     * @param A int[][] adjacency matrix
     */
    public void addPartition(int index, ArrayList<Integer> newPartition, int[][] A) {
        ArrayList<Integer> refinedPartition = new ArrayList<Integer>();
        if (allis1(newPartition)) {
            refinedPartition = newPartition;
        } else {
            refinedPartition = refinedPartitioning(newPartition, A[index]);
        }
        if (partitionList.size() == (index + 1)) {
            partitionList.add(refinedPartition);
        } else {
            partitionList.set(index + 1, refinedPartition);
        }
    }

    /**
     * Grund 3.3.8 In case if the row is canonical, with respect to former partition, we build
     * refined partitioning.
     *
     * @param partition ArrayList<Integer> atom partition
     * @param row int[] row
     * @return
     */
    public ArrayList<Integer> refinedPartitioning(ArrayList<Integer> partition, int[] row) {
        ArrayList<Integer> refined = new ArrayList<Integer>();
        int index = 0;
        int count = 1;
        for (Integer p : partition) {
            if (p != 1) {
                for (int i = index; i < p + index - 1; i++) {
                    if (i + 1 < p + index - 1) {
                        if (row[i] == row[i + 1]) {
                            count++;
                        } else {
                            refined.add(count);
                            count = 1;
                        }
                    } else {
                        if (row[i] == row[i + 1]) {
                            count++;
                            refined.add(count);
                            count = 1;
                        } else {
                            refined.add(count);
                            refined.add(1);
                            count = 1;
                        }
                    }
                }
                index = index + p;
            } else {
                index++;
                refined.add(1);
                count = 1;
            }
        }
        return refined;
    }

    public boolean blockTestPerm(
            int index,
            int r,
            int[][] A,
            ArrayList<Integer> partition,
            ArrayList<Integer> newPartition) {
        int y = findY(r);
        boolean check = true;
        int total = sum(partition);
        List<Permutation> cycleTrans = cycleTranspositions(index, partition);
        candidatePermutations(index, y, total, cycleTrans);
        boolean cycleCheck = check(index, y, total, A, partition, newPartition);
        if (!cycleCheck) {
            check = false;
        } else {
            addPartition(index, newPartition, A);
        }
        return check;
    }

    /** Build id permutation */
    public Permutation idPermutation(int size) {
        Permutation perm = new Permutation(size);
        return perm;
    }

    public int[] actArray(int[] strip, Permutation p) {
        int permLength = p.size();
        int arrayLength = strip.length;
        int[] modified = new int[arrayLength];
        for (int i = 0; i < permLength; i++) {
            modified[p.get(i)] = strip[i];
        }
        return modified;
    }

    public boolean equalBlockCheck(
            ArrayList<Integer> partition,
            int index,
            int y,
            int[][] A,
            Permutation cycleTransposition,
            Permutation perm) {
        boolean check = true;
        int[] canonical = A[index];
        int[] original = A[index];
        original =
                actArray(
                        actArray(A[cycleTransposition.get(index)], cycleTransposition),
                        perm); // siteden
        if (!Arrays.equals(canonical, original)) {
            check = false;
        }
        return check;
    }

    public int[] row2compare(int index, int[][] A, Permutation cycleTransposition) {
        int[] array = cloneArray(A[findIndex(index, cycleTransposition)]);
        array = actArray(array, cycleTransposition);
        return array;
    }

    public int findIndex(int index, Permutation cycle) {
        int size = cycle.size();
        int output = 0;
        for (int i = 0; i < size; i++) {
            if (cycle.get(i) == index) {
                output = i;
                break;
            }
        }
        return output;
    }

    public int[] cloneArray(int[] array) {
        int length = array.length;
        int[] cloned = new int[length];
        for (int i = 0; i < length; i++) {
            cloned[i] = array[i];
        }
        return cloned;
    }

    /**
     * Checking whether the modified and original rows are equal or not.
     *
     * @param cycleTransposition Permutation - cycle transposition
     * @param partition atom partition
     * @param index row index
     * @param y beginning index of a block
     * @param A adjacency matrix
     * @param testPermutation Permutation
     * @return
     */
    public int findIndex(int index, Permutation cycle, Permutation former) {
        int size = cycle.size();
        int output = 0;
        Permutation perm = cycle.multiply(former);
        for (int i = 0; i < size; i++) {
            if (perm.get(i) == index) {
                output = i;
                break;
            }
        }
        return output;
    }

    public int[] row2compare(
            int index, int[][] A, Permutation cycleTransposition, Permutation formerPermutation) {
        int[] array = cloneArray(A[findIndex(index, cycleTransposition, formerPermutation)]);
        Permutation perm = formerPermutation.multiply(cycleTransposition);
        array = actArray(array, perm);
        return array;
    }

    public Permutation getCanonicalPermutation(
            int[] max, int[] check, ArrayList<Integer> partition) {
        // size belirle
        int size = trialAbstractClass.size;
        int[] cycles = getCyclePermutation(max, check, partition);
        int[] perm = new int[size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i == cycles[j]) {
                    perm[i] = j;
                }
            }
        }
        return new Permutation(perm);
    }

    public int[] descendingSort(int[] array, int index0, int index1) {
        int temp = 0;
        for (int i = index0; i < index1; i++) {
            for (int j = i + 1; j < index1; j++) {
                if (array[i] < array[j]) {
                    temp = array[i];
                    array[i] = array[j];
                    array[j] = temp;
                }
            }
        }
        return array;
    }

    public int[] descendingSortWithPartition(int[] array, ArrayList<Integer> partition) {
        int i = 0;
        for (Integer p : partition) {
            array = descendingSort(array, i, i + p);
            i = i + p;
        }
        return array;
    }

    public boolean equalSetCheck(int[] original, int[] permuted, ArrayList<Integer> partition) {
        int[] temp = cloneArray(permuted);
        temp = descendingSortWithPartition(temp, partition);
        return equalSetCheck(partition, original, temp);
    }

    public Integer[] getBlocks(int[] row, int begin, int end) {
        return IntStream.range(begin, end).mapToObj(i -> row[i]).toArray(Integer[]::new);
    }

    public boolean equalSetCheck(ArrayList<Integer> partition, int[] original, int[] permuted) {
        boolean check = true;
        int i = 0;
        for (Integer p : partition) {
            Integer[] org = getBlocks(original, i, p + i);
            Integer[] perm = getBlocks(permuted, i, p + i);
            if (Arrays.equals(org, perm)) {
                i = i + p;
            } else {
                check = false;
                break;
            }
        }
        return check;
    }

    public int[] getCyclePermutation(int[] max, int[] check, ArrayList<Integer> partition) {
        int[] values = idValues(sum(partition));
        int i = 0;
        if (!equalSetCheck(max, check, partition)) {
            return values;
        } else {
            for (Integer p : partition) {
                Integer[] can = getBlocks(max, i, p + i);
                Integer[] non = getBlocks(check, i, p + i);
                values = getCyclesList(can, non, i, values);
                i = i + p;
            }
            return values;
        }
    }

    public int findMatch(Integer[] max, Integer[] non, int value, int start) {
        int size = non.length;
        int index = start;
        for (int i = start; i < size; i++) {
            if (non[i] == value) {
                if (max[i] != non[i]) {
                    index = i;
                    break;
                }
            }
        }
        return index;
    }

    public Integer[] permuteArray(Integer[] array, int i, int j) {
        int temp = 0;
        temp = array[i];
        array[i] = array[j];
        array[j] = temp;
        return array;
    }

    public int[] getCyclesList(Integer[] max, Integer[] non, int index, int[] values) {
        int i = 0;
        int permutationIndex = 0;
        while (i < max.length && max[i] != 0) {
            if (max[i] != non[i]) {
                permutationIndex = findMatch(max, non, max[i], i);
                if (i != permutationIndex) {
                    non = permuteArray(non, i, permutationIndex);
                }
                int temp = values[i + index];
                values[i + index] = values[permutationIndex + index];
                values[permutationIndex + index] = temp;
            }
            i++;
        }
        return values;
    }

    public Permutation getEqualPerm(
            Permutation cycleTransposition,
            int index,
            int y,
            int total,
            int[][] A,
            ArrayList<Integer> partition,
            ArrayList<Integer> newPartition) {
        int[] check = row2compare(index, A, cycleTransposition);
        Permutation canonicalPermutation = getCanonicalPermutation(A[index], check, newPartition);
        if (!biggerCheck(index, A[index], check, newPartition)) {
            biggest = false;
        }
        return canonicalPermutation;
    }

    public boolean biggerCheck(
            int index, int[] original, int[] permuted, ArrayList<Integer> partition) {
        permuted = descendingSortWithPartition(permuted, partition);
        // return descendingOrderUpperMatrix(index,partition, original, permuted);
        return descendingOrdercomparison(partition, original, permuted);
    }

    /**
     * Descending order check for a integer array
     *
     * @param array int[]
     * @return boolean
     */
    public boolean desOrderCheck(Integer[] array) {
        boolean check = true;
        int length = array.length;
        for (int i = 0; i < length - 1; i++) {
            if (array[i] < array[i + 1]) {
                check = false;
                break;
            }
        }
        return check;
    }

    /**
     * int array to int representation. This is needed for blockwise descending order check.
     *
     * @param array
     * @return
     */
    public int toInt(Integer[] array) {
        int result = 0;
        for (int i = 0; i < array.length; i++) {
            result = result * 10;
            result = result + array[i];
        }
        return result;
    }

    public boolean descendingOrdercomparison(
            ArrayList<Integer> partition, int[] original, int[] permuted) {
        boolean check = true;
        int i = 0;
        for (Integer p : partition) {
            Integer[] org = getBlocks(original, i, p + i);
            Integer[] perm = getBlocks(permuted, i, p + i);
            if (desOrderCheck(org)) {
                if (!Arrays.equals(org, perm)) {
                    if (toInt(org) == toInt(perm)) {
                        continue;
                    } else if (toInt(org) > toInt(perm)) {
                        check = true;
                        break; // If bigger than already in lexico order.
                    } else if (toInt(org) < toInt(perm)) {
                        check = false;
                        break;
                    }
                }
                i = i + p;
            } else {
                check = false;
                break;
            }
        }
        return check;
    }

    public Permutation getCanonicalCycle(
            int index,
            int y,
            int total,
            int[][] A,
            ArrayList<Integer> partition,
            ArrayList<Integer> newPartition,
            Permutation cycleTransposition) {
        equalAuto = true;
        biggest = true;
        Permutation canonicalPermutation = idPermutation(total);
        if (!equalBlockCheck(newPartition, index, y, A, cycleTransposition, canonicalPermutation)) {
            canonicalPermutation =
                    getEqualPerm(cycleTransposition, index, y, total, A, partition, newPartition);
            int[] check = row2compare(index, A, cycleTransposition);
            check = actArray(check, canonicalPermutation);
            if (canonicalPermutation.isIdentity()) {
                if (!descendingOrdercomparison(newPartition, A[index], check)) {
                    equalAuto = false;
                }
            }
        }
        return canonicalPermutation;
    }

    public boolean equalCheck(int[] former, int[] current, ArrayList<Integer> partition) {
        boolean check = true;
        int i = 0;
        for (int k = 0; k < partition.size(); k++) {
            Integer[] can = getBlocks(former, i, partition.get(k) + i);
            Integer[] org = getBlocks(current, i, partition.get(k) + i);
            if (!Arrays.equals(can, org)) {
                check = false;
                break;
            } else {
                i = i + partition.get(k);
                continue;
            }
        }
        return check;
    }

    public boolean check(
            int index,
            int y,
            int total,
            int[][] A,
            ArrayList<Integer> partition,
            ArrayList<Integer> newPartition) {
        boolean check = true;
        List<Permutation> formerList = new ArrayList<Permutation>();
        for (Permutation permutation : formerPermutations.get(index)) {
            Permutation canonicalPermutation =
                    getCanonicalCycle(index, y, total, A, partition, newPartition, permutation);
            if (biggest) {
                int[] test = row2compare(index, A, permutation);
                test = actArray(test, canonicalPermutation);
                if (descendingOrdercomparison(newPartition, A[index], test)) {
                    if (canonicalPermutation.isIdentity()) {
                        if (equalCheck(A[index], test, partition)) {
                            formerList.add(permutation);
                        }
                    } else {
                        Permutation newPermutation = canonicalPermutation.multiply(permutation);
                        formerList.add(newPermutation);
                    }
                } else {
                    formerList.clear();
                    check = false;
                    break;
                }
            } else {
                formerList.clear();
                check = false;
                break;
            }
        }

        if (check) {
            formerPermutations.get(index).clear();
            formerPermutations.set(index, formerList);
        }
        return check;
    }

    public void candidatePermutations(int index, int y, int total, List<Permutation> cycles) {
        List<Permutation> newList = new ArrayList<Permutation>();
        newList.addAll(cycles);
        formerPermutations.add(index, newList);
        if (index != 0) {
            List<Permutation> formers = formerPermutations.get(index - 1);
            for (Permutation form : formers) {
                if (!form.isIdentity()) {
                    formerPermutations.get(index).add(form);
                }
            }
            List<Permutation> newForm = new ArrayList<Permutation>();
            for (Permutation frm : formers) {
                if (!frm.isIdentity()) {
                    newForm.add(frm);
                }
            }
            List<Permutation> newCycles = new ArrayList<Permutation>();
            for (Permutation cyc : cycles) {
                if (!cyc.isIdentity()) {
                    newCycles.add(cyc);
                }
            }
            for (Permutation perm : newForm) {
                for (Permutation cycle : newCycles) {
                    Permutation newPermutation = cycle.multiply(perm);
                    if (!newPermutation.isIdentity()) {
                        formerPermutations.get(index).add(newPermutation);
                    }
                }
            }
        }
    }

    public int sum(ArrayList<Integer> list, int index) {
        int sum = 0;
        for (int i = 0; i <= index; i++) {
            sum = sum + list.get(i);
        }
        return sum;
    }

    /**
     * Summing entries of a list.
     *
     * @param list List<Integer>
     * @return int sum
     */
    public int sum(List<Integer> list) {
        int sum = 0;
        for (int i = 0; i < list.size(); i++) {
            sum = sum + list.get(i);
        }
        return sum;
    }

    public void clearFormers(boolean check, int y) {
        if (check == false) {
            for (int i = y; i < formerPermutations.size(); i++) {
                formerPermutations.get(i).removeAll(formerPermutations.get(i));
            }
            for (int i = y + 1; i < partitionList.size(); i++) {
                partitionList.get(i).removeAll(partitionList.get(i));
            }
        }
    }

    public boolean blockTest(int r, int[][] A) {
        boolean check = true;
        int y = findY(r);
        int z = findZ(r);
        for (int i = y; i <= z; i++) {
            if (!blockTestPerm(
                    i, r, A, partitionList.get(i), canonicalPartition(i, partitionList.get(i)))) {
                check = false;
                break;
            }
        }
        clearFormers(check, y);
        return check;
    }
}
