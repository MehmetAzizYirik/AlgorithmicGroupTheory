/**
 * MIT License
 *
 * Copyright (c) 2018 Mehmet Aziz Yirik
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


/**
 * This class is the re-implementation of group algorithms from C.A.G.E.S. book.[1] 
 * The functions are from Chapter 6: Groups and Symmetry. The source codes, given with 
 * the book, have been converted from C to Java.
 * 
 * [1] Kreher, D.L. and Stinson, D.R., 1998. Combinatorial algorithms: generation,
 * enumeration, and search (Vol. 7). CRC press.
 * 
 * @author Mehmet Aziz Yirik
 */


package Orbits;

import java.util.ArrayList;
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
