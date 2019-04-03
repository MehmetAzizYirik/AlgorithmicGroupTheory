package Orbits;

import java.util.Arrays;
import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
import org.openscience.cdk.group.PermutationGroup.Backtracker;


public class norb {
		public static PermutationGroup group=PermutationGroup.makeSymN(4);
		public static boolean DoneEarly;
		public static int X;
		public static int[] T5= new int[group.getSize()+1];
		public static int[] T6= new int[group.getSize()+1];
		public static int[] T7= new int[group.getSize()+1];
		public static boolean[] T8 = new boolean[group.getSize()+1];
		
		public static int choose(int n, int k) { 
			int C[][] = new int[n+1][k+1]; 
			int i, j;
			for (i = 0; i <= n; i++){ 
				for (j = 0; j <= Math.min(i, k); j++){ 
					if (j == 0 || j == i) {
						C[i][j] = 1; 
					}else {
						C[i][j] = C[i-1][j-1] + C[i-1][j];
					}
				}
			} 
			return C[n][k]; 
	    } 
		
		public static void type(Permutation g, int[] T) {
			int i,j,ell;
			int size=group.getSize();
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
			if (backtracker.isFinished()) { 
	            return;
	        }
	        if (l == group.getSize()) {
	            backtracker.applyTo(g);
	        } else {
	            for (int i = 0; i < group.getSize(); i++) {
	                Permutation h = group.get(l,i);
	                if (h != null) { 
	                	backtrack(l + 1, g.multiply(h), backtracker);
	                }
	            }
	        }
	    }
		
		public static void nOrb(PermutationGroup G, int[] N, int limit) {
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
		
		public static void main(String[] args) {
			PermutationGroup groupp= PermutationGroup.makeSymN(4);
			int[] N= new int[group.getSize()+1];
			nOrb(groupp,N,4);
		}
	}
