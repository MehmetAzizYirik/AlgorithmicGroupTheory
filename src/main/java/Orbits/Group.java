package Orbits;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;

public class Group {
	 public static boolean DoneEarly=false;
	 public static int i;
	 public static int X=0;
	 public static int size=4;	 
	 public static Permutation I =  newPerm(); 
	 public static Permutation T0 = newPerm();
	 public static Permutation T1 = newPerm();
	 public static Permutation T2 = newPerm();
	 public static Permutation T3 = newPerm();
	 public static Permutation C1 = newPerm();
	 public static Permutation C2 = newPerm();
	 public static int[] T5 = new int[size+1];
	 public static int[] T6 = new int[size+1];
	 public static int[] T7 = new int[size+1];
	 public static boolean[] T8 = new boolean[size+1];
	 public static List<Permutation> T4= newPermList(size);  
	 public static ArrayList<Permutation> getPi(){
		 ArrayList<Permutation> list=new ArrayList<Permutation>();
		 for(i=0;i<(size+1);i++) {
			 list.add(newPerm());
		 }
		 return list;
	 }
	 public static ArrayList<Permutation> pi= getPi();
	 public static Permutation newPerm() {
		 Permutation p = new Permutation(size);
		 return p;
	 } 
	 public static PermutationGroup buildGroup(int n) {
		 PermutationGroup group=PermutationGroup.makeSymN(n);
		 return group;
	 }
	 public static List<Permutation> newPermList(int m){
		 List<Permutation> L= new ArrayList<Permutation>(m);
		 for(int k=0;k<m;k++) {
			 L.add(null);
		 }
		 return(L);
	 }
}
