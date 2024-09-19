/******************************************************************************************************************************************
 * If you use this code, please cite
 * L. Cirigliano, C. Castellano and G. Bianconi
  * "General theory for extended-range percolation on simple and multiplex networks
  * Phys. Rev. E 110, 034302 (2024)
 *  ***************************************************************************************************************************************
 *This code implements Extended Range Percolation (ERP) on single Poisson networks with N nodes and average degree c
 *The ERP has range N and the Monte Carlo simulations are averaged over Nrunmax iterations.range
 *
 *
 *INPUT
 N total number of nodes
 c average degree of Poisson layers
 R Range of extended range percolation
 Nrunmax average over which the MERGC is calculated.

 
 OUTPUT
 
 The output file "ERP_Monte_Carlo_simulations.txt" contains three columns:
 p=1-f the fraction of untruested nodes  and the order patermeters of  ERP
 given by  P_infty, and U_infty obtained by direct numerical calculations
 * */




#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 20000
#define c  3.0
#define R  4
#define Nrunmax 100
 

#define err_tol 0.01
int *vis1,*size_cluster1,**knn1,*k1,c1,c2,c3,*occ,*dam1,*drange;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component in each layer
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster,int rm){
	int j, n3,aus_cluster_size,n,j2,aus;

    if(dam1[i]==1) cluster_size++;
        vis1[i]=ncluster;
        for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
                if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster,1);
				cluster_size=aus_cluster_size;
                }
            if((dam1[j]==0)&&(drange[j]>rm)&&(rm<R)){
                drange[j]=rm;
                    aus_cluster_size=Recurrence(j, cluster_size, ncluster,rm+1);
                    cluster_size=aus_cluster_size;
                }
		}
	
		return cluster_size;
}

int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GC,nn,mu,*occ2,cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc,np,GCU;
	int s1,s2,Nc1,Nc2,Nc3,aus,aus3,**adj1,**adj2,**adj3,N0,i_aus;
	float p,x,f,**xd1,*Sm,**n1,MGC1,nsum1,nsumold1,nsum2,aus1;
    int *GCm,*GCUm;
	
	   char filename[60],string1[50],string2[50];
	
	FILE *gp3,*gp;
    


  /**************************************************************************
   open file for output 
   finemane GC
   at the end of teh program the file will contain 
   two columns: p GC GC2
   
   **************************************************************************/

		srand48(time(NULL));
	

	N0=N;
	

	vis1=(int*)calloc(N,sizeof(int));
	occ=(int*)calloc(N,sizeof(int));
	k1=(int*)calloc(N,sizeof(int));
	xd1=(float**)calloc(N,sizeof(float*));
	dam1=(int*)calloc(N,sizeof(int));
    drange=(int*)calloc(N,sizeof(int));
	knn1=(int**)calloc(N,sizeof(int*));
    n1=(float**)calloc(N,sizeof(float*));
    GCm=(int*)calloc(100,sizeof(int));
    GCUm=(int*)calloc(100,sizeof(int));
		for(i=0;i<N;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
            n1[i]=(float*)calloc(N,sizeof(float));
        }
    for(i=0;i<N;i++){
			xd1[i]=(float*)calloc(N,sizeof(float));
		}
	size_cluster1=(int*)calloc(N,sizeof(int));
	
    
    for(i=0;i<N;i++){
        k1[i]=0;
    }

	for(i=0;i<N;i++){
        for (j=i+1;j<N;j++){
            x=drand48();
            if(x<c/(float)N){
                k1[i]++;
                k1[j]++;
            
                knn1[i][k1[i]-1]=j;
                knn1[j][k1[j]-1]=i;
            
                n1[i][j]=(int)(2*drand48());
                n1[j][i]=(int)(2*drand48());
                nsum1+=n1[i][j]+n1[j][i];
        }
    }
    }
		
  

				
/*%%%%%%%%%%%%SIMULATION%%%%%%%%%%%%%%%%%*/
	for (nrun=0;nrun<Nrunmax;nrun++){	
	
        for(i=0;i<N;i++){
            xd1[i][0]=drand48();
            dam1[i]=1;
            drange[i]=R+1;
            
        }
    
        for(nc=0;nc<100;nc++){
		f=nc*0.01;
            
            for(i=0;i<N;i++){
                dam1[i]=1;
                drange[i]=R+1;
                
            }
           
                for (i=0;i<N;i++){
                 if(xd1[i][0]<f){
                     dam1[i]=0;}
                 else {dam1[i]=1;}
             }
            
		for(i=0;i<N;i++){
			vis1[i]=0;
		}
    
		GC=0;
        m1=0;
		ncluster1=0;
		for(n=0;n<N;n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(n, cluster_size,ncluster1,1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}				
			}			
        }
        Nc1=c1;
        GC=0;
        GCU=0;
        for(i=0;i<N;i++){
            if((vis1[i]==Nc1)&&(dam1[i]==1)){
                GC++;
				//	m2++;
            }
        }
            for(i=0;i<N;i++){
                if((vis1[i]==Nc1)&&(dam1[i]==0)){
                    GCU++;
                    //    m2++;
                }
            }
        GCm[nc]+=GC;
            GCUm[nc]+=GCU;
        }
    }
    gp3=fopen("ERP_Monte_Carlo_simulations.txt","w");
        for (nc=0;nc<100;nc++){
            f=nc*0.01;
        
        fprintf(gp3,"%f %f %f\n",1-f,(float)GCm[nc]/(float)(Nrunmax*N),(float)GCUm[nc]/(float)(Nrunmax*N));
        }
		fclose(gp3);
	
	return 0;	
}

