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
 * The functions are from Chapter 6: Groups and Symmetry. Page 216-222. The source codes,
 * given with the book, have been converted from C to Java.
 * 
 * [1] Kreher, D.L. and Stinson, D.R., 1998. Combinatorial algorithms: generation,
 * enumeration, and search (Vol. 7). CRC press.
 * 
 * @author Mehmet Aziz Yirik
 */



package Orbits;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.openscience.cdk.group.Permutation;
import org.openscience.cdk.group.PermutationGroup;
//TODO: The error is I used the only one V value array but it should be like given in the class form. like a->u->V.
public class orbitRepresentation {
	public static PermutationGroup group=PermutationGroup.makeSymN(4);
	public static int size= group.getSize();
	public static boolean DoneEarly;
	public static int[] V= new int[size];
	public static int WORDSIZE = 32;
	public static int FIRSTBIT(int x) {
		return (31-((x) & 037777600000 ? ((x) & 037700000000 ? \ leftbit[((x)>>24) & 0377] : 8+leftbit[(x)>>16]) \ : ((x) & 0177400 ? 16+leftbit[(x)>>8] : 24+leftbit[x])))
	}
	//public static int LASTBIT(x) wherebit[(x&(-x))%37]

	
	public static Set<Integer> GetSet(Set<Integer> A){
		Set<Integer> B=  new HashSet<Integer>(A);
		return B;
	}
	
	public static Set GetEmptySet(){ //What do they mean with getemptyset ?
		Set emptySet = Collections.emptySet();
		int i;
		for(i=0;i<WORDS;i++) V[i]=0;
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
	
	public static int FindLargestElement(Set S){
		int i,j,k,a;
		i=group.getSize();
		while(!V[i]) i--;
		a=V[i];
		j=LASTBIT(a);
		k=(WORDSIZE-1-j)+(i*WORDSIZE); // Not clear why they need WORDSIZE ?
		return( (j>=0)? k : -1 );
	}
	public static int FindPrevElement(int u, Set S) {
	/*
	**   Returns the largest element in the set S
	*/
		int i,j,k;
		int a;

		j=(u % WORDSIZE );
		i=u/WORDSIZE;
		a=V[i]&~rightmask[j];

		if(!a){
			i--;
			while(!V[i]) i--; 
			a=V[i];
		}
		j=LASTBIT(a);
		k=(WORDSIZE-1-j)+(i*WORDSIZE);
		return( (j>=0)? k : -1 );
	}
	
	public static int SetOrder(Set S){
		/*
		**  Algorithm 1.8
		**
		** Return the number of ones in the set with
		** SubsetLexRank S.
		*/
		int ans,i;
		int x;
		ans = 0;
		for(i=0;i<WORDS;i++){
		  x=V[i];
		  while (x){
			  ans += look[x&(UINT)0377];
			  x >>= (UINT)8;
		  }
		}
		return(ans);
	}
	public static void Apply(Permutation g,Set<Integer> A,Set<Integer> B){
		int u,h,z;
		Set<Integer> S2 = new HashSet<Integer>();
		S2=GetEmptySet(); // Describe S2
		h=0;
		u=FindLargestElement(A);
		z=SetOrder(A);
		while(h<z){
			h=h+1;
			SetInsert(V[u],S2);
			if(h<z) u=FindPrevElement(u,A);
		}
		GetSet(B);
	}
	
	public static void MinRepBT(int k,PermutationGroup G,int ell){
		int i,j,x,m,r,start,n;
		//setlist C;
		Permutation g=null;
		n=G.getSize();
		if (ell==0){
			start=0;
		}else{
			start=X[ell-1];
		}
		m=n;
		for(x=start;x<n;x++) { 
			if(MemberOfSet(x,blocks[ell])){
				r=0;
				while((r<ell)&&(X[r]==OptX[r])) r=r+1;
				if( (r<ell) && X[r]>=OptX[r]) return;
				X[ell]=x;
				SetPerm(ell+1,X,beta);
				ChangeBase(G,beta);
				for(j=0;j<=x;j++) {
					if((g==T[ell][j])!=null) break;
				}
				if (V[x]<=m){
					m=V[x];
					X[ell]=m;
					if (k==ell+1){
						for(i=r;(i<k) && (X[i]==OptX[i]);i++) {
							if((i!=k) && ORB->X[i]<ORB->OptX[i]) { 
								for(i=0;i<k;i++) { 
									OptX[i]=X[i];
								}
							}
						}
					}
					else{
						Apply(g,blocks[ell],blocks[ell+1]);
						SetDelete(m,blocks[ell+1]);
						MinRepBT(k,G,ell+1);
					}
				}
			}
		}
	}
	
	
	public static void MinRep(PermutationGroup G, Set A, Set ans) {
	/*
	**  Algorithm 6.15
	*/
		int i,k,deg;
		//setlist C;
		deg=G.getSize();
		k=0; 
		for(i=0;i<deg;i++){
			if(MemberOfSet(i,A)) {
				OptX[k++]=i;
			}
		}
		GetSet(blocks[0],A);
		MinRepBT(k,G,0);
		ans=GetEmptySet(ans);
		for(i=0;i<k;i++) {
			SetInsert(OptX[i],ans);
		}
	}
	
	public static void main(String[] args) {
		

	}

}
