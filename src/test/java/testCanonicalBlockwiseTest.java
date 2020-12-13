import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.interfaces.IAtomContainer;

import AlgorithmicGroupTheory.canonicalBlockwiseTest;

public class testCanonicalBlockwiseTest extends canonicalBlockwiseTest{
	
	/**
	 * The matrix is from Grund Thesis - Example 3.3.12, 3.6.4
	 */
	
	public static int[][] mat= new int[10][10];
	public static void initializeMatrix() {
		mat[0][2]=1;
		mat[0][3]=1;
		mat[0][4]=1;
		mat[1][2]=1;
		mat[1][3]=1;
		mat[1][5]=1;
		mat[2][0]=1;
		mat[2][1]=1;
		mat[2][4]=2;
		mat[3][0]=1;
		mat[3][1]=1;
		mat[3][5]=1;
		mat[3][6]=1;
		mat[4][0]=1;
		mat[4][2]=2;
		mat[4][7]=1;
		mat[5][1]=1;
		mat[5][3]=1;
		mat[5][8]=1;
		mat[5][9]=1;
		mat[6][3]=1;
		mat[7][4]=1;
		mat[8][5]=1;
		mat[9][5]=1;
	}
	
	/**
	 * Initializing an atom partition from Grund Example 3.3.12
	 */
	
	public static ArrayList<Integer> partition= new ArrayList<Integer>();
	public static void initializePartition() {
		partition.add(2);
		partition.add(4);
		partition.add(4);
	}
	
	public static int[][] hydrogenFreeMatrix= new int[10][10];
	public static void initializehydrogenFreeMatrix() {
		hydrogenFreeMatrix[0][2]=1;
		hydrogenFreeMatrix[0][3]=1;
		hydrogenFreeMatrix[0][4]=1;
		hydrogenFreeMatrix[1][2]=1;
		hydrogenFreeMatrix[1][3]=1;
		hydrogenFreeMatrix[1][5]=1;
		hydrogenFreeMatrix[2][0]=1;
		hydrogenFreeMatrix[2][1]=1;
		hydrogenFreeMatrix[2][4]=2;
		hydrogenFreeMatrix[3][0]=1;
		hydrogenFreeMatrix[3][1]=1;
		hydrogenFreeMatrix[3][5]=1;
		hydrogenFreeMatrix[4][0]=1;
		hydrogenFreeMatrix[4][2]=2;
		hydrogenFreeMatrix[5][1]=1;
		hydrogenFreeMatrix[5][3]=1;
	}
	
	/**
	 * Grund Example 3.3.12. The hydrogen columns start from the 6th one.
	 */
	
	@Test
	public void testAddHydrogens() {
		initializeMatrix();
		initializehydrogenFreeMatrix();
		
		ArrayList<Integer> degree= new ArrayList<Integer>();
		degree.add(3);
		degree.add(3);
		degree.add(4);
		degree.add(4);
		degree.add(4);
		degree.add(4);
		
		initialDegrees=degree;
		Assert.assertArrayEquals(addHydrogens(hydrogenFreeMatrix, 6),mat);
	}
	
	@Test
	public void testAddRepresentatives() {
		addRepresentatives(0,new Permutation(1,2,0,3,4));
		Assert.assertEquals(representatives.size(),1);
	}
	
	/**
	 * Building an atom container from connectivity matrix
	 * @throws CloneNotSupportedException 
	 */
	
	@Test 
	public void testBuildC() throws CloneNotSupportedException {
		initializeMatrix();
		IAtomContainer atomContainer = buildC(mat);
		Assert.assertEquals(atomContainer.getBondCount(),12);
	}
	
	@Test
	public void testCanonicalBlockbasedGenerator() {
		//TODO: How to test recursive function ?
	}
	
	@Test
	public void testCanonicalPartition() {
		initializePartition();
		ArrayList<Integer> canonicalPartition= new ArrayList<Integer>();
		canonicalPartition.add(1);
		canonicalPartition.add(1);
		canonicalPartition.add(4);
		canonicalPartition.add(4);
		Assert.assertEquals(canonicalPartition(0, partition),canonicalPartition);
	}
	
	@Test
	public void testCInverse() {
		//TODO: I need to read it again.I dont know what it outputs.
	}
	
	/**
	 * Clearing representatives and partitions of a block.
	 */
	
	@Test
	public void testClear() {
		initializePartition();
		partitionList.add(partition);
		addRepresentatives(0,new Permutation(0,1,2,3,4));
		addRepresentatives(0,new Permutation(1,2,0,4,3));
		/**
		 * remove the entries of the first block from representatives
		 * and partitionList.
		 */
		clear(true,0);
		Assert.assertEquals(representatives.size(),0);
		Assert.assertEquals(partitionList.size(),0);
	}
	
	@Test
	public void testClone() throws CloneNotSupportedException, IOException, CDKException {
		initializeMatrix();
		atomContainer=buildC(mat);
		int atoms= atomContainer.getAtomCount();
		Assert.assertEquals(clone(atomContainer).getAtomCount(),atoms);
	}
	
	/**
	 * Summing the entries in a column until the given final index.
	 */
	@Test
	public void testCsum() {
		initializeMatrix();
		Assert.assertEquals(Csum(mat, 0, 1, 6),3);
	}
	
	/**
	 * Summing the entries in a row until the given final index.
	 */
	@Test
	public void testLsum() {
		initializeMatrix();
		Assert.assertEquals(Lsum(mat, 0, 0, 6),3);
	}
	
	@Test
	public void testCycleTransposition() {
		/**
		 * Grund Example 3.3.4
		 */
		ArrayList<Integer> partition= new ArrayList<Integer>();
		partition.add(2);
		partition.add(1);
		partition.add(2);
	
		Assert.assertEquals(cycleTranspositions(0, partition),2);
	}
	
	@Test
	public void testDegreeSeq() {
		/**
		 * Get the degree sequence from a molecular formula.
		 */
		ArrayList<Integer> degrees= new ArrayList<Integer>();
		degrees.add(4);
		degrees.add(4);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);

		Assert.assertEquals(degreeSeq("C2H6"),degrees);
	}
	
	@Test
	public void testDepict() {
		//TODO: How to test depiction function ?
	}
	
	@Test
	public void testDescBlockCheck() {
		initializeMatrix();
		/**
		 * Grund 3.3.19
		 */
		ArrayList<Integer> partition = new ArrayList<Integer>();
		partition.add(1);
		partition.add(1);
		partition.add(1);
		partition.add(1);
		partition.add(1);
		partition.add(1);
		partition.add(4);
		boolean check= descBlockCheck(new Permutation(0,1,3,2,4,5,6,7,8,9), partition, 2, 1, mat, new Permutation(0,1,2,3,4,5,6,7,8,9));
		Assert.assertEquals(check, true);
	}
	
	/**
	 * Without permutation action, checking whether the row
	 * is in the descendng form or not.
	 */
	
	@Test 
	public void testDescBlockwiseCheck2() {
		initializeMatrix();
		initializePartition();
		
		/**
		 * Grund Example 3.3.12: Testing the first row with its canonical
		 * partition
		 */
		
		boolean check= descBlockwiseCheck2(0, 0, mat, canonicalPartition(0, partition));
		Assert.assertEquals(check, true);
	}
	
	/**
	 * Descending order check for a Integer array 
	 */
	
	@Test
	public void testDesOrderCheck() {
		Integer[] array= new Integer[3];
		array[0]=1;
		array[1]=0;
		array[2]=1;
		boolean check=desOrderCheck(array);
		Assert.assertEquals(check, false);
	}
	
	@Test
	public void testDistributeHydrogens() throws FileNotFoundException, UnsupportedEncodingException, CloneNotSupportedException, CDKException {
		/**
		 * Distribute hydrogens to the carbons in C6H6.
		 */
		
		ArrayList<Integer> partition = new ArrayList<Integer>();
		partition.add(6);
		partition.add(6);
		
		ArrayList<Integer> degrees = new ArrayList<Integer>();
		degrees.add(4);
		degrees.add(4);
		degrees.add(4);
		degrees.add(4);
		degrees.add(4);
		degrees.add(4);
		
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		degrees.add(1);
		
		List<ArrayList<Integer>> hydrogens= distributeHydrogens(partition, degrees);
		Assert.assertEquals(hydrogens,7);
	}
	
	/**
	 * Building an atom partition of a molecular formula.
	 */
	@Test
	public void testDistribution() {
		ArrayList<Integer> partition = new ArrayList<Integer>();
		partition.add(6);
		partition.add(6);
		
		Assert.assertEquals(distribution("C6H6"),partition);
	}
	
	/**
	 * Checks whether the rows are identical after permutation actions.
	 */
	
	@Test
	public void testEqualBlockCheck() {
		/**
		 * Grund Example 3.3.19, first row with
		 * cycle (0,1) and canonical permutation (4,5).
		 */
		
		initializeMatrix();
		ArrayList<Integer> partition = new ArrayList<Integer>();
		partition.add(1);
		partition.add(1);
		partition.add(4);
		partition.add(4);
		
		boolean check= equalBlockCheck(new Permutation(1,0,2,3,4,5,6,7,8,9),partition, 0,0 ,mat, new Permutation(0,1,2,3,5,4,6,7,8,9));
		Assert.assertEquals(check,true);
	}
	
	/**
	 * X value is always 0 for the first block since there is no former
	 * block in the system.
	 */
	
	@Test
	public void testFindX() {
		initializePartition();
		inputPartition=partition;
		//For the first block of the partition
		Assert.assertEquals(findX(0),0);
	}
	
	/**
	 * Y value of the first block is always 0.
	 */
	
	@Test
	public void testFindY() {
		initializePartition();
		inputPartition=partition;
		//For the first block of the partition
		Assert.assertEquals(findY(0),0);
	}
	
	/**
	 * Z value is the last row index in a block.
	 */
	
	@Test
	public void testFindZ() {
		initializePartition();
		inputPartition=partition;
		//For the first block of the partition
		Assert.assertEquals(findZ(0),1);
	}
	
	@Test
	public void testFormerBlocksRepresentatives() {
		//TODO: This function will be changed, so no need of a test.  
	}
	
	/**
	 * Bundan sonra sadece basit olanlari kodladim; digerleri bir ornek uzerinden 
	 * yapilmali.
	 */
	
	/**
	 * Getting a part of a integer array
	 */
	
	@Test
	public void testGetBlock() {
		int[] array = new int[5];
		array[0]=1;
		array[1]=2;
		array[2]=3;
		array[3]=4;
		array[4]=5;
		
		Integer[] output = new Integer[4];
		output[0]=1;
		output[1]=2;
		output[2]=3;
		output[3]=4;
		
		Assert.assertArrayEquals(getBlocks(array,0,3),output);
		
	}
	
	@Test
	public void testPartitionCriteria() {
		initializePartition();
		ArrayList<Integer> output = new ArrayList<Integer>();
		output.add(1);
		output.add(1);
		output.add(4);
		output.add(4);
		
		/**
		 * The first re-partitioning 
		 */
		
		Assert.assertEquals(partitionCriteria(partition,1),output);
	}
	
	/**
	 * In an 5*5 matrix, getting the former entry's indices.
	 */
	
	@Test
	public void testPredecessor() {
		ArrayList<Integer> indices = new ArrayList<Integer>();
		indices.add(1);
		indices.add(4);
		
		ArrayList<Integer> formerIndices = new ArrayList<Integer>();
		formerIndices.add(1);
		formerIndices.add(3);
		
		Assert.assertEquals(predecessor(indices,5), formerIndices);
	}
	
	/**
	 * In an 5*5 matrix, getting the next entry's indices.
	 */
	
	@Test
	public void testSuccessor() {
		ArrayList<Integer> indices = new ArrayList<Integer>();
		indices.add(1);
		indices.add(3);
		
		ArrayList<Integer> nextIndices = new ArrayList<Integer>();
		nextIndices.add(1);
		nextIndices.add(4);
		
		Assert.assertEquals(successor(indices,5), nextIndices);
	}
	
	/**
	 * After hydrogen distribution, removing the hydrogens from the partition 
	 */
	
	@Test
	public void testRedefineThePartition() {
		initializePartition();
		
		ArrayList<Integer> output = new ArrayList<Integer>();
		output.add(2);
		output.add(4);
		
		Assert.assertEquals(redefineThePartition(partition), output);
	}
	
	@Test
	public void testRefinedPartitioning() {
		initializePartition();
		
		int[] array = new int[10];
		array[0] = 1;
		array[1] = 0;
		array[2] = 1;
		array[3] = 1;
		array[4] = 0;
		array[5] = 0;
		array[6] = 0;
		array[7] = 0;
		array[8] = 0;
		array[9] = 0;
		
		/**
		 * Refining a partition with respect to the entries 
		 * in a row.
		 */
		
		ArrayList<Integer> output = new ArrayList<Integer>();
		output.add(1);
		output.add(1);
		output.add(2);
		output.add(2);
		output.add(4);
		
		Assert.assertEquals(refinedPartitioning(partition, array),output);
	}
	
	/**
	 * Getting the entry indices based on the atom partition
	 */
	
	@Test
	public void testSubPartitionsList() {
		initializePartition();
		
		ArrayList<Integer> first = new ArrayList<Integer>();
		first.add(0);
		first.add(1);
		
		ArrayList<Integer> second = new ArrayList<Integer>();
		second.add(2);
		second.add(3);
		second.add(4);
		second.add(5);
		
		ArrayList<Integer> third = new ArrayList<Integer>();
		third.add(6);
		third.add(7);
		third.add(8);
		third.add(9);
		
		ArrayList<ArrayList<Integer>> output = new ArrayList<ArrayList<Integer>>();
		
		output.add(first);
		output.add(second);
		output.add(third);
		
		Assert.assertEquals(subPartitionsList(partition),output);
	}
	
	@Test
	public void testSum() {
		int[] array = new int[3];
		array[0]=1;
		array[1]=1;
		array[2]=1;
		
		Assert.assertEquals(sum(array),3);			
	}
	
	@Test 
	public void testToInt() {
		Integer[] array = new Integer[3];
		array[0]=1;
		array[1]=0;
		array[2]=2;
		
		Assert.assertEquals(toInt(array),102);
	}
	
	/**
	 * In the backwardBuild function, we might need to backward
	 * to the former block, then r needs to be updated. We check
	 * indices position in the initial partition.  
	 */
	
	@Test
	public void testUpdateR() {
		ArrayList<Integer> indices = new ArrayList<Integer>();
		indices.add(1);
		indices.add(4);
		
		/**
		 * The first index is smaller than the Y value ,2, of the current 
		 * block. Then we are already moved to the former block. r=0.
		 */
		
		Assert.assertEquals(updateR(1,indices),0);
	}
	
	/**
	 * Checks all entries are zero or not.
	 */
	
	@Test
	public void testAllIs0(){
		ArrayList<Integer> list= new ArrayList<Integer>();
		list.add(1);
		list.add(0);
		list.add(0);
		Assert.assertEquals(list,false);
	}
	
	@Test
	public void testBackwardBuild() {
		//TODO: How can we test a recursive function ?
	}
	
	@Test
	public void testBlock() {
		//TODO: ?
	}
	
	@Test
	public void testBlockTest() {
		//TODO: ?
	}
	
	/**
	 * Constructing an atom container
	 * @throws CDKException 
	 * @throws CloneNotSupportedException 
	 * @throws IOException 
	 */
	@Test
	public void testBuild() throws IOException, CloneNotSupportedException, CDKException {
		build("C6H6");
		Assert.assertEquals(atomContainer.getAtomCount(),12);
	}
	
	
	/**
	 * General Functions
	 */
	
	@Test
	public void testaddOnes() {
		int number=5;
		ArrayList<Integer> list= new ArrayList<Integer>();
		addOnes(list,number);
		Assert.assertEquals(list.size(),number);
	}
	
	/**
	 * Grund - Example 3.3.12
	 */
	
	@Test
	public void testaddPartition() {
		initializeMatrix();
		ArrayList<Integer> partition = new ArrayList<Integer>();
		partition.add(1);
		partition.add(1);
		partition.add(4);
		partition.add(4);
		addPartition(0, partition, mat);
		Assert.assertEquals(partitionList.size(),1);
	}
	
	
	/**
	 * Connectivity Test
	 * 
	 * Testing the functions from Grund's Thesis - 3.6.2
	 */
	
	@Test
	public void testnValues() {
		initializeMatrix();
		Set<Integer> nSet= new HashSet<Integer>();
		nSet.add(0);
		nSet.add(2);
		nSet.add(3);
		nSet.add(4);
		
		/**
		 * Testing the first row for the nValues.
		 * From the row 6, the rows are for atoms with degree 1.
		 */
		
		Assert.assertEquals(nValues(0, 10, mat),nSet);
	}
	
	@Test
	public void testwValues() {
		ArrayList<Integer> kFormer= new ArrayList<Integer>();
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(6);
		kFormer.add(7);
		kFormer.add(8);
		kFormer.add(9);
		
		Set<Integer> nSet= new HashSet<Integer>();
		nSet.add(3);
		nSet.add(5);
		nSet.add(6);
		
		Set<Integer> wSet = new HashSet<Integer>();
		wSet.add(0);
		wSet.add(6);
		
		Assert.assertEquals(wValues(nSet, kFormer),wSet);
	}
	
	@Test
	public void testkValues() {
		ArrayList<Integer> kFormer= new ArrayList<Integer>();
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(0);
		kFormer.add(8);
		kFormer.add(9);
		
		Set<Integer> wSet = new HashSet<Integer>();
		wSet.add(0);
		wSet.add(8);
		wSet.add(9);
		
		ArrayList<Integer> output= new ArrayList<Integer>();
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		output.add(0);
		
		Assert.assertEquals(kValues(10, wSet,kFormer), output);
	}
		
	@Test
	public void testConnectivityTest() {
		int[][] matrix = new int[10][10];
		matrix[0][1]=1;
		matrix[0][2]=2;
		matrix[1][0]=1;
		matrix[1][2]=2;
		matrix[2][0]=2;
		matrix[2][1]=2;
		matrix[3][4]=3;
		matrix[3][5]=1;
		matrix[4][3]=3;
		matrix[4][6]=1;
		matrix[5][3]=1;
		matrix[5][7]=1;
		matrix[5][8]=1;
		matrix[5][9]=1;
		matrix[6][4]=1;
		matrix[7][5]=1;
		matrix[8][5]=1;
		matrix[9][5]=1;
		
		/**
		 * The adjacency matrix represents a molecules with 2 
		 * saturated sub components. 
		 */
		
		Assert.assertEquals(connectivityTest(6, matrix),false);
	}

}
