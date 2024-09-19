/******************************************************************************************************************************************
 * If you use this code, please cite
 * L. Cirigliano, C. Castellano and G. Bianconi
  * "General theory for extended-range percolation on simple and multiplex networks"
  * Phys. Rev. E 110, 034302 (2024)
 *  ***************************************************************************************************************************************
 * This program generates a  multiplex network with M=2 Poisson layers with the same average degree c
 * t consider a Nrunmax random configuration of initial damage
 * and evaluates the  Extended Range Mutually Connected Giant Component as a function of the fraction of  nodes damaged
 
 INPUT
 N total number of nodes in each layer
 c average degree of Poisson layers
 R Range of extended range percolation
 Version: Option 1 (version A) and option 2 (version B)
 Nrunmax average over which the MERGC is calculated.

 
 OUTPUT
 
 The output file "ERP_MCGC_MonteCarlo.txt"  contains three columns:
 p=1-f the fraction of untruested nodes  and the order patermeters of  MERP
 given by  P_infty, and U_infty obtained by direct numerical calculations

 
*************************************************************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#define N  20000
#define c  2.2
#define R  2
#define Version 2
#define Nrunmax  50


int *vis1,*vis2,*vis3,*size_cluster1,*size_cluster2,*size_cluster3,**knn1,**knn2,**knn3,*k1,*k2,*k3,c1,c2,c3,*occ,*dam1,*dam2,*dam3,*drange;


/***********************************************************************************************
 Recurrence is a subrutine for calculating the giant component in each layer
 *************************************************************************************************/

int Recurrence( int net,int i , int cluster_size, int ncluster,int rm){
    int j, n3,aus_cluster_size,n,j2,aus;
    if(net==1){
    if(dam1[i]==1) cluster_size++;
        vis1[i]=ncluster;
        for(n3=0;n3<k1[i];n3++){
            j=knn1[i][n3];
                if((vis1[j]==0)&&(dam1[j]==1)){
                aus_cluster_size=Recurrence(1,j, cluster_size, ncluster,1);
                cluster_size=aus_cluster_size;
                }
            if((dam1[j]==0)&&(drange[j]>rm)&&(rm<R)){
                drange[j]=rm;
                    aus_cluster_size=Recurrence(1,j, cluster_size, ncluster,rm+1);
                    cluster_size=aus_cluster_size;
                }
        }
    }
    if(net==2){
    if(dam2[i]==1) cluster_size++;
        vis2[i]=ncluster;
        for(n3=0;n3<k2[i];n3++){
            j=knn2[i][n3];
                if((vis2[j]==0)&&(dam2[j]==1)){
                aus_cluster_size=Recurrence(2,j, cluster_size, ncluster,1);
                cluster_size=aus_cluster_size;
                }
            if((dam2[j]==0)&&(drange[j]>rm)&&(rm<R)){
                drange[j]=rm;
                    aus_cluster_size=Recurrence(2,j, cluster_size, ncluster,rm+1);
                    cluster_size=aus_cluster_size;
                }
        }
    }
        return cluster_size;
}

/*
int Recurrence(int net, int i , int cluster_size, int ncluster, int rm){
	int j, n3,aus_cluster_size;
	
	cluster_size++;	
	if(net==1){
		vis1[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(1,j, cluster_size, ncluster,1);
				cluster_size=aus_cluster_size;
			}
	
		}
	}
	if(net==2){
		vis2[i]=ncluster;
		for(n3=0;n3<k2[i];n3++){
			j=knn2[i][n3];			
			if((vis2[j]==0)&&(dam2[j]==1)){
			aus_cluster_size=Recurrence(2,j, cluster_size, ncluster,1);
			cluster_size=aus_cluster_size;
			}
		}
	}
	

	return cluster_size;
}*/
/*******************************************************************************************
 RecurrenceM is a subrutine for calulating the  Mutually Connected Giant Component
 *******************************************************************************************/

int RecurrenceM( int net, int i , int cluster_size, int ncluster,int rm){
    int j, n3,aus_cluster_size,n,j2,aus;
    if(net==1){
    if(dam1[i]==1) cluster_size++;
        occ[i]=ncluster;
        for(n3=0;n3<k1[i];n3++){
            j=knn1[i][n3];
                if((occ[j]==0)&&(dam1[j]==1)&&(vis2[j]==c2)&&(dam2[j]==1)){
                aus_cluster_size=RecurrenceM(1,j, cluster_size, ncluster,1);
                cluster_size=aus_cluster_size;
                }
            if((dam1[j]==0)&&(drange[j]>rm)&&(rm<R)){
                if((Version==2)||((Version==1)&&(vis2[j]==c2))){
                    drange[j]=rm;
                    aus_cluster_size=RecurrenceM(1,j, cluster_size, ncluster,rm+1);
                    cluster_size=aus_cluster_size;
                    }
            }
        }
    }
    if(net==2){
    if(dam2[i]==1) cluster_size++;
        occ[i]=ncluster;
        for(n3=0;n3<k2[i];n3++){
            j=knn2[i][n3];
                if((occ[j]==0)&&(dam2[j]==1)&&(vis1[j]==c1)&&(dam1[j]==1)){
                aus_cluster_size=RecurrenceM(2,j, cluster_size, ncluster,1);
                cluster_size=aus_cluster_size;
                }
            if((dam2[j]==0)&&(drange[j]>rm)&&(rm<R)){
                if((Version==2)||((Version==1)&&(vis1[j]==c1))){
                    drange[j]=rm;
                    aus_cluster_size=RecurrenceM(2,j, cluster_size, ncluster,rm+1);
                    cluster_size=aus_cluster_size;
                }
            }
        }
    }
                   
        return cluster_size;
}


/* RecurrenceM(int net, int i , int cluster_size, int ncluster, int rm){
	int j, n3,aus_cluster_size;
	float *xd1,*xd2,*xd3;
	
	cluster_size++;
	
	if(net==1){
		occ[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if(dam1[j]==1){
			if(((vis2[j]==c2)&&(occ[j]==0)&&(dam2[j]==1))){
				aus_cluster_size=RecurrenceM(1,j, cluster_size, ncluster,1);
				cluster_size=aus_cluster_size;
			}
			}
		}
	}
	if(net==2){
		occ[i]=ncluster;
		for(n3=0;n3<k2[i];n3++){
			j=knn2[i][n3];		
			if(dam2[j]==1){
			if(((vis1[j]==c1)&&(occ[j]==0)&&(dam1[j]==1))){
				aus_cluster_size=RecurrenceM(2,j, cluster_size, ncluster,1);
				cluster_size=aus_cluster_size;
			}
		}
		}
	}
	


	return cluster_size;
}*/



int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GMCC, cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc;
	int s1,s2,Nc1,Nc2,Nc3,**n1,**n2,**n3,aus,aus1,aus2,aus3,**adj1,**adj2,**adj3,RMCGC,nsum,nsumold,layer,maxdenR,**denR,n_aus,MCGCU;
	int *occ1old,*occ2old,*occ3old,*S11_1,*S11_2,*S11_3,*S1_1,*S1_2,*S1_3,S11,MCGC;
	float p,x,f,*xd1,*xd2,*xd3,*Sm,RMCGCm,S11m,*MCGCav,*MCGCUav;

	
    char filec[60],**s,string1[50],string2[50],filename[60];
	
	FILE *gp2,*gp,*ffp,*gp4;
  

  /**************************************************************************
   open file for output 
  **************************************************************************/

		srand48(time(NULL));

	
	Sm=(float*)calloc(200,sizeof(float));
	vis1=(int*)calloc(N,sizeof(int));
	vis2=(int*)calloc(N,sizeof(int));
	vis3=(int*)calloc(N,sizeof(int));
	occ1old=(int*)calloc(N,sizeof(int));
	occ2old=(int*)calloc(N,sizeof(int));
	occ3old=(int*)calloc(N,sizeof(int));
	S11_1=(int*)calloc(N,sizeof(int));
	S11_2=(int*)calloc(N,sizeof(int));
	S11_3=(int*)calloc(N,sizeof(int));
	S1_1=(int*)calloc(N,sizeof(int));
	S1_2=(int*)calloc(N,sizeof(int));
	S1_3=(int*)calloc(N,sizeof(int));
	occ=(int*)calloc(N,sizeof(int));
	k1=(int*)calloc(N,sizeof(int));
	k2=(int*)calloc(N,sizeof(int));
	k3=(int*)calloc(N,sizeof(int));
	xd1=(float*)calloc(N,sizeof(float));
	xd2=(float*)calloc(N,sizeof(float));
	xd3=(float*)calloc(N,sizeof(float));
	dam1=(int*)calloc(N,sizeof(int));
	dam2=(int*)calloc(N,sizeof(int));
	dam3=(int*)calloc(N,sizeof(int));
	sigma=(int*)calloc(N,sizeof(int));
	knn1=(int**)calloc(N,sizeof(int*));
	knn2=(int**)calloc(N,sizeof(int*));
	knn3=(int**)calloc(N,sizeof(int*));
	n1=(int**)calloc(N,sizeof(int*));
	n2=(int**)calloc(N,sizeof(int*));
	n3=(int**)calloc(N,sizeof(int*));
	adj1=(int**)calloc(N,sizeof(int*));
	adj2=(int**)calloc(N,sizeof(int*));
	adj3=(int**)calloc(N,sizeof(int*));
    MCGCav=(float*)calloc(100,sizeof(float));
    MCGCUav=(float*)calloc(100,sizeof(float));
	denR=(int**)calloc(101,sizeof(int*));
    drange=(int*)calloc(N,sizeof(int));
	for(i=0;i<101;i++){
		denR[i]=(int*)calloc(N,sizeof(int));
	}
		for(i=0;i<N;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
			knn2[i]=(int*)calloc(N,sizeof(int));
			knn3[i]=(int*)calloc(N,sizeof(int));
			adj1[i]=(int*)calloc(N,sizeof(int));
			adj2[i]=(int*)calloc(N,sizeof(int));
			adj3[i]=(int*)calloc(N,sizeof(int));
			n1[i]=(int*)calloc(N,sizeof(int));
			n2[i]=(int*)calloc(N,sizeof(int));
			n3[i]=(int*)calloc(N,sizeof(int));

		}
	size_cluster1=(int*)calloc(N,sizeof(int));
	size_cluster2=(int*)calloc(N,sizeof(int));
	size_cluster3=(int*)calloc(N,sizeof(int));

	
	for(i=0;i<N;i++){
		k1[i]=0;
		k2[i]=0;
	}
	nsum=0;
	
	/***************************************************************************
	 Generate  a  multiplex network
	 **************************************************************************/
	for (i=0;i<N;i++){
			for (j=i+1;j<N;j++){
					  if(drand48()<(c/((float)N))){
					  k1[i]++;
					  k1[j]++;
					  knn1[i][k1[i]-1]=j;
					  knn1[j][k1[j]-1]=i;
					  
					  adj1[i][j]=1;
					  adj1[j][i]=1;
				
					  }
					  
					  if(drand48()<(c/((float)N))){
					  k2[i]++;
					  k2[j]++;
					  knn2[i][k2[i]-1]=j;
					  knn2[j][k2[j]-1]=i;
					  
					  adj2[i][j]=1;
					  adj2[j][i]=1;
				
					  }
			}
	}
	
					   		
		
	
	for(nrun=0;nrun<Nrunmax;nrun++){
	
	for(i=0;i<N;i++){
		xd1[i]=drand48();	}
	

		for(nc=0;nc<100;nc++){
			f=nc*0.01;
			p=1-f;
		for(i=0;i<N;i++){
			if(xd1[i]<f){
				dam1[i]=0;
                dam2[i]=0;
            }
            else {dam1[i]=1;dam2[i]=1;}
			
        }
				
		
		
		
		/***************************************************************************************
		 Calculates the giant component of each layer
		 ****************************************************************************************/
		ncluster1=0;
		ncluster2=0;
		
		for(i=0;i<N;i++){
			vis1[i]=0;
			vis2[i]=0;
			
		}
		m1=0;
		m2=0;
		
            for(i=0;i<N;i++){

                drange[i]=R+1;

            }
		for(n=0;n<N;n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(1,n, cluster_size, ncluster1,1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}
				
			}
        }
            for(i=0;i<N;i++){

                drange[i]=R+1;

            }
            for(n=0;n<N;n++){
			if((vis2[n]==0)&&(dam2[n]==1)){
				cluster_size=0;
				ncluster2++;
				cluster_size=Recurrence(2,n, cluster_size, ncluster2,1);
				size_cluster2[ncluster2]=cluster_size;
				if (cluster_size>m2) {m2=cluster_size;c2=ncluster2;}
				
			}
			
		
		}
		
		
		
		Nc2=ncluster2;
		Nc1=ncluster1;
		
		/*************************************************************************************
		 Calculates back and forth the giant components in each layer until the RMCGC is found 
		 GMCC is the size of the RMCGC
		 *************************************************************************************/
		GMCC=0;
		m1_aus=0;
		m2_aus=0;
		n_aus=0;
		//m3_aus=0;
		while(m1_aus!=m1){
			n_aus++;
			m1=m1_aus;
			m2=m2_aus;
			m1_aus=0;
			m2_aus=0;

			for(i=0;i<N;i++)occ[i]=0;
            for(i=0;i<N;i++){
                drange[i]=R+1;
            }
			for(i=0;i<N;i++){
				if(dam2[i]==1){
					if(((vis1[i]==c1)&&(vis2[i]==c2)&&(occ[i]==0)&&(dam1[i]==1))){
						cluster_size=0;
						ncluster2++;
						cluster_size=RecurrenceM(2,i, cluster_size, ncluster2,1);
						if (cluster_size>m2_aus) {m2_aus=cluster_size;c2_aus=ncluster2;}
					}				
				}
			}
			c2=c2_aus;

			for(i=0;i<N;i++){vis2[i]=occ[i];occ[i]=0;}
            for(i=0;i<N;i++){
                drange[i]=R+1;
            }
			
			n_aus++;
			for(i=0;i<N;i++){
				if(dam1[i]==1){
					if(((vis1[i]==c1)&(vis2[i]==c2)&(occ[i]==0)&&(dam2[i]==1))){
						cluster_size=0;
						ncluster1++;
						cluster_size=RecurrenceM(1,i, cluster_size, ncluster1,1);
						if (cluster_size>m1_aus) {m1_aus=cluster_size;c1_aus=ncluster1;}
					}				
				}
			}
			c1=c1_aus;
			for(i=0;i<N;i++){vis1[i]=occ[i];occ[i]=0;}

		}
            MCGCU=0;
            for(i=0;i<N;i++){
                if((dam1[i]==0)&&((vis1[i]==c1)||(vis2[i]==c2))){
                    MCGCU++;}
            }
            
		MCGC=m1;
        MCGCav[nc]+=((float)m1/(float)N)/(float)Nrunmax;
        MCGCUav[nc]+=((float)MCGCU/(float)N)/(float)Nrunmax;
        }
	}
    sprintf(filename,"ERP_MCGC_MonteCarlo.txt");
    gp4=fopen(filename,"w");
	for (nc=0;nc<100;nc++){
		p=1-nc*0.01;
		fprintf(gp4,"%lf %lf %lf \n",p,MCGCav[nc],MCGCUav[nc]);
    }
    fclose(gp4);
	
	
	

			
	return 0;	
}

