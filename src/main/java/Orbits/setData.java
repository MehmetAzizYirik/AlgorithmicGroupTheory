package Orbits;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.openscience.cdk.group.Permutation;

public class setData {
	public static int n=4;
	ArrayList<Integer> FullSet = newSet();
	ArrayList<Integer> EmptySet= newSet();
	 
	public static ArrayList<Integer> newSet() {
		ArrayList<Integer> list=new ArrayList<Integer>(n);
		return list;
	}
}
