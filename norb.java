package Orbits;

/**
This is the re-implementation of the norb (number of orbits) function given in C.A.G.E.S book by Kreher.
**/

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.group.PermutationGroup.Backtracker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class norb {
		public static PermutationGroup group=PermutationGroup.makeSymN(4);
		//public static int WORDS= computeWords(group.getSize());
		public static int WORDSIZE=32;
		public static boolean DoneEarly;
		public static int X;
		public static void setPerm(int[] V) {
			Permutation perm= new Permutation(group.getSize());
			V=perm.getValues();
		}
		public static int[] T5= new int[group.getSize()+1];
		public static Permutation perm=new Permutation(group.getSize());
		//public static int[] T5= perm.getValues();
		class group_data{
			Permutation I;
			Permutation T0;
			Permutation T1;
			Permutation T2;
			Permutation T3;
			Permutation C1;
			Permutation C2;
			Permutation T4;
			int T6;
			int T7;
			boolean[] T8 = new boolean[100];
			Permutation pi;
		}
		
		class orb_data{
			Set S;
			Set S2;
			setlist C= new setlist();
			int X;
			int OptX;
			Permutation beta;
			Permutation f;
			Permutation Invf;
			Permutation F;
		}
		
		class set_data{
			Set FullSet;
			Set EmptySet;
		}
		
		class setlist{
			universe U = new universe();
			int size ;
			Set blocks;
		}
		
		class iso_data{
			setlist TQ = new setlist();
			setlist TP0 = new setlist(); 
			setlist TP1 = new setlist();
			setlist TP2 = new setlist();
			setlist TP3 = new setlist();
			setlist TR = new setlist();
			Permutation best;
			int BestExist;
			Set Tx;
			Set Ty;
			setlist L= new setlist();
			Set U;
			setlist S= new setlist();
			Set T;
			int N;
			Permutation f;
			Permutation Invf;
		}
		
		class universe{
			int NODES;
			int WORDS;
			int n=group.getSize();
			group_data GRP= new group_data();
			orb_data ORB = new orb_data(); 
			set_data SET = new set_data(); 
			iso_data ISO = new iso_data();
		}
		public static int computeWords(int n) {
			return (n/WORDSIZE+1);
		}
		
		public static int choose2(int n, int k) {
			int k2, ans;
			if(n < k) {
				return 0;
			}
			k2 = n - k < k ? n - k : k;
			ans = 1;
			for(int i = 1; i <= k2; i++) {
				ans = (ans * (n + 1 - i)) / i;
			}
			return ans;
		}
		
		static int choose(int n, int k) 
	    { 
	    int C[][] = new int[n+1][k+1]; 
	    int i, j; 
	      
	        // Calculate  value of Binomial Coefficient in bottom up manner 
	    for (i = 0; i <= n; i++) 
	    { 
	        for (j = 0; j <= Math.min(i, k); j++) 
	        { 
	            // Base Cases 
	            if (j == 0 || j == i) 
	                C[i][j] = 1; 
	       
	            // Calculate value using previosly stored values 
	            else
	                C[i][j] = C[i-1][j-1] + C[i-1][j]; 
	          } 
	     } 
	       
	    return C[n][k]; 
	    } 
		
		public static group_data groupData;
		public static IChemObjectBuilder builder =SilentChemObjectBuilder.getInstance();
		public static boolean[] T8 = new boolean[group.getSize()+1];
		public static void type(Permutation g, int[] T) {
			int i,j,ell;
			//Universe of the group is same for all the permutations of the group.
			//Describe the universe in the group first. That is the same for all the group members. You dont have different ones.
			int size=group.getSize();
			//boolean[] p= new boolean[size];		
			for(i=0;i<size;i++) {
				T8[i]=true;
				T[i+1]=0;
			}
			for(i=0;i<size;i++) {
				if(T8[i]) {
					ell=1;
					j=i;
					T8[j]=false;
					while(T8[g.getValues()[j]]) {
						ell=ell+1;
						j=g.getValues()[j];
						T8[j]=false;
					}
					T[ell]=T[ell]+1;
				}
			}
		}
		public static int[] T6= new int[group.getSize()+1];
		public static int[] T7= new int[group.getSize()+1];
		public static void RecPartition(int k,int m,int B,int N){ 
			int[] c= new int[1000];
			int i,prod,n;
			n=group.getSize();
			if(m==0){
				for(i=1;i<=n;i++) {
					c[i]=0;
				}
				for(i=1;i<=N;i++) {
					c[T6[i]]=c[T6[i]]+1;
				}
				prod=1; 
				for(i=1;i<=n;i++) {
					prod=prod*choose(T7[i],c[i]);
				}
				T5[k]=T5[k]+prod;
			}else{
				for(i=1;i<=Math.min(B,m);i++){
					T6[N+1]=i;
					RecPartition(k,m-i,i,N+1);
				}
			}
		}
		
		public static void norb(PermutationGroup G, int[] N, int limit) {
			int k,deg,OrderG;
			deg=group.getSize();
			DoneEarly=false;
			X=limit;
			for(k=0;k<=deg;k++) T5[k]=0;
			norbBack();
			OrderG=Math.toIntExact(group.order()); 
			System.out.println(Arrays.toString(T5));
			for(k=0;k<=deg;k++) N[k]=T5[k]/OrderG;
			System.out.println(Arrays.toString(N));
		}
		
		
		public static void norbBack() {
	        Backtracker norbUse = new Backtracker() {
	        	public void applyTo (Permutation g) {
	        		type(g,T7);
	        		for(int k=0;k<=X;k++){
	        			T6[1]=k;
	        			RecPartition(k,k,k,0);
	        		}
	        	}
	            public boolean isFinished() {
	                return false;
	            }
	        };
	        apply(norbUse);
	    }
		
		public static void apply(Backtracker backtracker) {
	        backtrack(0, new Permutation(group.getSize()), backtracker);
	    }
		
		public static void backtrack(int l, Permutation g, Backtracker backtracker) {
			if (backtracker.isFinished()) { // He coded this backtrack in general. If needed finish is checked.
	            return;
	        }
	        if (l == group.getSize()) {
	            backtracker.applyTo(g);
	        } else {
	            for (int i = 0; i < group.getSize(); i++) {
	                Permutation h = group.get(l,i);
	                if (h != null) { // Because he defined them as null. Need to understand sims chain
	                	//In the book it says mult but the mult is done in backtrack function g.multiply(h)
	                	backtrack(l + 1, g.multiply(h), backtracker);
	                }
	            }
	        }
	    }
		
		public static Set<Integer> GetSet(Set<Integer> A){
			Set<Integer> B=  new HashSet<Integer>(A);
			return B;
		}
		
		public static Set GetEmptySet(){
			Set emptySet = Collections.emptySet();
			return emptySet;
		}
		
		public static void SetInsert(int u, Set<Integer> A){
			A.add(u);
		}
		
		public static boolean MemberOfSet(int u, Set<Integer> A) {
			boolean check = false;
			if(A.contains(u)) {
				check=true;
			}
			return check;
		}
		public static void SetPerm(int k,int[] A,Permutation g){
			int i,j,n=group.getSize();
			Set empty= GetEmptySet();
			int[] values=g.getValues();
			for(i=0;i<k;i++){
				SetInsert(A[i],empty);
				values[i]=A[i];
			}
			for(j=0;j<n;j++) if(!MemberOfSet(j,empty)){
				values[i]=j;
				i=i+1;
			}
		}
		
		public static void main(String[] args) {
			// TODO Auto-generated method stub
			PermutationGroup groupp= PermutationGroup.makeSymN(4);
			int[] N= new int[group.getSize()+1];
			norb(groupp,N,4);
		}
	}
}
