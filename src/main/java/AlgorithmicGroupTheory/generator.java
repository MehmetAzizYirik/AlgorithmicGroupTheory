package AlgorithmicGroupTheory;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.io.SDFWriter;

public class generator {
    public int size = 0;
    public int hIndex = 0;
    public int count = 0;
    public int matrixSize = 0;
    public SDFWriter outFile;
    public String filedir;
    public boolean flag = true;
    public ArrayList<int[][]> output = new ArrayList<int[][]>();
    public PrintWriter pWriter;
    public boolean stripIterate = true;
    public boolean entryChanges = false;
    public int[][] max;
    public int[][] L;
    public int[][] C;
    public int[] degrees;

    generator(int[] degree2) {
        this.degrees = degree2;
    }

    public Map<String, Integer> valences;

    {
        // The atom valences are from CDK.
        valences = new HashMap<String, Integer>();
        valences.put("C", 4);
        valences.put("N", 3);
        valences.put("O", 2);
        valences.put("S", 2);
        valences.put("P", 3);
        valences.put("F", 1);
        valences.put("I", 1);
        valences.put("Cl", 1);
        valences.put("Br", 1);
        valences.put("H", 1);
    }
    /** Basic functions */

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

    /**
     * Summing entries of an array.
     *
     * @param array int[]
     * @return int sum
     */
    public int sum(int[] array) {
        int sum = 0;
        for (int i = 0; i < array.length; i++) {
            sum = sum + array[i];
        }
        return sum;
    }

    public int sum(ArrayList<Integer> list, int index) {
        int sum = 0;
        for (int i = 0; i <= index; i++) {
            sum = sum + list.get(i);
        }
        return sum;
    }

    /** Candidate Matrix Generation Functions */

    /**
     * L; upper triangular matrix like given in 3.2.1. For (i,j), after the index, giving the
     * maximum line capacity.
     *
     * @param degrees
     * @return upper triangular matrix
     */
    public void upperTriangularL() {
        L = new int[hIndex][hIndex];
        for (int i = 0; i < hIndex; i++) {
            for (int j = i + 1; j < hIndex; j++) {
                L[i][j] = Math.min(degrees[i], Lsum(i, j + 1));
            }
        }
    }

    /**
     * C; upper triangular matrix like given in 3.2.1. For (i,j), after the index, giving the
     * maximum column capacity.
     *
     * @param degrees
     * @return upper triangular matrix
     */
    public int[][] upperTriangularC() {
        C = new int[hIndex][hIndex];
        for (int i = 0; i < hIndex; i++) {
            for (int j = i + 1; j < hIndex; j++) {
                C[i][j] = Math.min(degrees[j], Csum(i + 1, j));
            }
        }
        return C;
    }

    public int Lsum(int i, int j) {
        int sum = 0;
        for (int k = j; k < hIndex; k++) {
            sum = sum + max[i][k];
        }
        return sum;
    }

    public int Csum(int i, int j) {
        int sum = 0;
        for (int k = i; k < hIndex; k++) {
            sum = sum + max[k][j];
        }
        return sum;
    }

    /**
     * Theorem 3.2.1
     *
     * @param degrees degree sequence
     */
    public void maximalMatrix() {
        max = new int[hIndex][hIndex];
        for (int i = 0; i < hIndex; i++) {
            for (int j = 0; j < hIndex; j++) {
                int di = degrees[i];
                int dj = degrees[j];
                if (i == j) {
                    max[i][j] = 0;
                } else {
                    if (di != dj) {
                        max[i][j] = Math.min(di, dj);
                    } else if (di == dj && i != j) {
                        max[i][j] = (di - 1);
                    }
                }
            }
        }
    }

    /**
     * public void genStrip(int[] degreeList, int y, int z) throws IOException,
     * CloneNotSupportedException { int[][] A = new int[matrixSize][matrixSize]; degrees=degreeList;
     * flag=true; maximalMatrix(); upperTriangularL(); upperTriangularC(); int[] indices= new
     * int[2]; indices[0]=0; indices[1]=1; boolean callForward=true;
     *
     * <p>/** r y z needs to be calculated from the initialPartition before the canonical test.
     *
     * <p>but how will we know where we reach the block end ? so we need these values in the
     * extension step not something to calculate before canonicalTest.
     */

    /**
     * while(flag) { int[][] mat=nextStep(A,indices,callForward); if(!flag) { break; }
     * indices=successor(indices,max.length); callForward=true; nextStep(mat,indices,callForward); }
     * }*
     */
    public int[] successor(int[] indices, int size) {
        int i0 = indices[0];
        int i1 = indices[1];
        if (i1 < (size - 1)) {
            indices[0] = i0;
            indices[1] = (i1 + 1);
        } else if (i0 < (size - 2) && i1 == (size - 1)) {
            indices[0] = (i0 + 1);
            indices[1] = (i0 + 2);
        }
        return indices;
    }

    public int[] predecessor(int[] indices, int size) {
        int i0 = indices[0];
        int i1 = indices[1];
        if (i0 == i1 - 1) {
            indices[0] = (i0 - 1);
            indices[1] = (size - 1);
        } else {
            indices[0] = i0;
            indices[1] = (i1 - 1);
        }
        return indices;
    }

    /**
     * public int[][] nextStep(int[][] A, int[] indices,boolean callForward) throws IOException,
     * CloneNotSupportedException{ if(callForward) { return forwardRow(A, indices,callForward);
     * }else { return backward(A,indices,callForward); } }*
     */

    /**
     * After generating matrices, adding the hydrogen with respect to the pre-hydrogen distribution.
     *
     * @param A int[][] adjacency matrix
     * @param index int beginning index for the hydrogen setting
     * @return
     */
    public int[][] addHydrogens(int[][] A, int index) {
        int hIndex = index;
        int hydrogen = 0;
        for (int i = 0; i < index; i++) {
            hydrogen = trialAbstractClass.initialDegrees[i] - sum(A[i]);
            for (int j = hIndex; j < (hIndex + hydrogen); j++) {
                A[i][j] = 1;
                A[j][i] = 1;
            }
            if (hydrogen != 0) {
                hIndex = hIndex + hydrogen;
            }
        }
        return A;
    }

    /**
     * In backward function, updating the block index.
     *
     * @param r int block index
     * @param indices ArrayList<Integer> atom valences
     * @return
     */
    public void updateR(int[] indices, int y, int z, int r) {
        if (indices[0] < y) {
            r--;
        } else if (indices[0] > z) {
            r++;
        }
    }

    /**
     * The third line of the backward method in Grund 3.2.3. The criteria to decide which function
     * is needed: forward or backward.
     *
     * @param x the value in the adjacency matrix A[i][j]
     * @param lInverse lInverse value of indices {i,j}
     * @param l
     * @return
     */
    public boolean backwardCriteria(int x, int lInverse, int l) {
        boolean check = false;
        int newX = (x - 1);
        if ((lInverse - newX) <= l) {
            check = true;
        }
        return check;
    }

    /**
     * The third step in Grund 3.2.3.
     *
     * <p>Backward step in the algorithm.
     *
     * @param A int[][] adjacency matrix
     * @param indices ArrayList<Integer> indices
     * @throws IOException
     * @throws CloneNotSupportedException
     * @throws CDKException
     */

    /**
     * public int[][] backward(int[][] A, int[] indices, boolean callForward) throws IOException,
     * CloneNotSupportedException { int i=indices[0]; int j=indices[1]; if(i==0 && j==1) {
     * flag=false; return A; }else { indices=predecessor(indices, max.length); //updateR(indices, y,
     * z, r); i= indices[0]; j= indices[1]; int x= A[i][j]; int l2= LInverse(degrees,i,j,A); int c2=
     * CInverse(degrees,i,j,A); /** I changed in backcriteria from x to x-1 but then I had error
     * c6h6 was 98 not 217
     */
    /**
     * if(x>0 && (backwardCriteria((x),l2,L[i][j]) && backwardCriteria((x),c2,C[i][j]))){
     * A[i][j]=(x-1); A[j][i]=(x-1); indices = successor(indices,max.length); //updateR(indices);
     * callForward=true; return nextStep(A,indices,callForward); }else { callForward=false; return
     * nextStep(A,indices,callForward); } } }*
     */

    /**
     * public int[][] forwardRow(int[][] A, int[] indices, boolean callForward) throws IOException,
     * CloneNotSupportedException { int i=indices[0]; int j=indices[1]; int lInverse=
     * LInverse(degrees,i,j,A); int cInverse= CInverse(degrees,i,j,A); int minimal =
     * Math.min(max[i][j],Math.min(lInverse,cInverse));
     *
     * <p>/** First step in the forward method.
     */

    /**
     * int maximumValue = forwardMaximal(minimal, lInverse, L[i][j], cInverse, C[i][j]);
     * callForward=true; if(j==(max.length-1)) { int[][] newMat= forwardSubRow(lInverse, cInverse,
     * maximumValue, i, j,A,indices,callForward); return newMat; }else { return
     * forwardSubRow(lInverse, cInverse, maximumValue, i, j,A,indices,callForward); } }*
     */

    /**
     * public int[][] forwardSubRow(int lInverse, int cInverse, int maximalX,int i, int j, int[][]
     * A, int[] indices, boolean callForward) throws CloneNotSupportedException, IOException {
     * if(((lInverse-maximalX)<=L[i][j]) && ((cInverse-maximalX)<=C[i][j])) { A[i][j]=maximalX;
     * A[j][i]=maximalX; if(i==(max.length-2) && j==(max.length-1)) { if(canonicalTest(A)) {
     * //if(isCanonical(A)) { output.add(A); //A=addHydrogens(A,hIndex); int[][] mat2= new
     * int[A.length][A.length]; for(int k=0;k<A.length;k++) { for(int l=0;l<A.length;l++) {
     * mat2[k][l]=A[k][l]; } } //IAtomContainer molden= buildC(addHydrogens(mat2,hIndex));
     * if(connectivityTest(hIndex,addHydrogens(mat2,hIndex))){ System.out.println("canonical
     * matrix"+" "+count); for(int s=0;s<6;s++) { for(int k=0; k<s;k++) { System.out.print(" "); }
     * for(int k=s+1; k<6;k++) { System.out.print(mat2[s][k]); } System.out.println(); }
     *
     * <p>/**pWriter.println("canonical matrix"+" "+count); for(int s=0;s<6;s++) { for(int k=0;
     * k<s;k++) { pWriter.print(" "); } for(int k=s+1; k<6;k++) { pWriter.print(mat2[s][k]); }
     * pWriter.println(); }*
     */
    // IAtomContainer mol= buildC(addHydrogens(mat2,hIndex));
    // outFile.write(mol);
    // depict(mol,filedir+count+".png");
    /**
     * count++; } } callForward=false; return nextStep(A, indices, callForward); }else {
     * if(indices[0]==findZ(r) && indices[1]==(max.length-1)) { callForward=canonicalTest(A);
     * //callForward=isCanonical(A); if(callForward) { indices=successor(indices,max.length);
     * updateR(indices); } return nextStep(A, indices, callForward); }else { if(i<findZ(r) &&
     * j==(max.length-1)) { return A; }else { indices=successor(indices,max.length);
     * //updateR(indices); return nextStep(A, indices, true); } } } }else { callForward=false;
     * return nextStep(A, indices,callForward); } }*
     */

    /**
     * First step in the forward method.
     *
     * @param min minimum of L, C amd maximal matrices for {i,j} indices.
     * @param lInverse Linverse value of {i,j}
     * @param l L value of {i,j}
     * @param cInverse Cinverse value of {i,j}
     * @param c C value of {i,j}
     * @return int max
     */
    public int forwardMaximal(int min, int lInverse, int l, int cInverse, int c) {
        int max = 0;
        for (int v = min; v >= 0; v--) {
            if (((lInverse - v) <= l) && ((cInverse - v) <= c)) {
                max = max + v;
                break;
            }
        }
        return max;
    }

    /**
     * Grund Thesis 3.2.1 - Page 35 L and C Inverse formulas given at the bottom of the page.
     *
     * @param degrees atom valences
     * @param i row index
     * @param j column index
     * @param A matrix in the generation process
     * @return
     */
    public int LInverse(int[] degrees, int i, int j, int[][] A) {
        int sum = 0;
        for (int s = 0; s < j; s++) {
            sum = sum + A[i][s];
        }
        return degrees[i] - sum;
    }

    public int CInverse(int[] degrees, int i, int j, int[][] A) {
        int sum = 0;
        for (int s = 0; s < i; s++) {
            sum = sum + A[s][j];
        }
        return degrees[j] - sum;
    }
}
