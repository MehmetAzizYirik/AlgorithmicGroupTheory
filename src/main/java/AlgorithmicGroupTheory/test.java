package AlgorithmicGroupTheory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;

import AlgorithmicGroupTheory.PermutationGroupFunctions;

public class test extends PermutationGroupFunctions {
	/**
	 * Grund Thesis, example 3.6.4.
	 */
	public static int[][] mat= new int[10][10];
	public static void initialMat() {
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
	 * Connectivity Test
	 * 
	 * Testing the functions from Grund's Thesis - 3.6.2
	 */
	
	@Test
	public void testnValues() {
		initialMat();
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
	public void testconnectivityTest() {
		initialMat();
		Assert.assertEquals(connectivityTest(6, mat),true);
	}
	
	public static void main(String[] args) throws CloneNotSupportedException, CDKException, IOException { 
		
		mat[0][1]=1;
		mat[0][2]=2;
		mat[1][0]=1;
		mat[1][2]=2;
		mat[2][0]=2;
		mat[2][1]=2;
		mat[3][4]=3;
		mat[3][5]=1;
		mat[4][3]=3;
		mat[4][6]=1;
		mat[5][3]=1;
		mat[5][7]=1;
		mat[5][8]=1;
		mat[5][9]=1;
		mat[6][4]=1;
		mat[7][5]=1;
		mat[8][5]=1;
		mat[9][5]=1;
		System.out.println(connectivityTest(6, mat));
	}
}
