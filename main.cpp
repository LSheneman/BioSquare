//
//  main.cpp
//  BioSquareGT
//
//  Original code Created by Arend Hintze on 2/6/13, modified and extended by Jory Schossau.
//  Copyright (c) Arend Hintze All rights reserved.
//


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <map>
#include <math.h>
#include <string.h>
#include "memoryUsageLib.h"

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif


#define randDouble ((double)rand()/(double)RAND_MAX)

const int GENES=2;//10;

int globalUpdate=0;
using namespace std;
double w=1.0; //selection strength
double mutationRate=0.01;
double replacementRate=0.01;
int generations=500000; // these are updates, but we call them generations which is a misnomer
bool brh=false;

int xm[8]={0,1,1,1,0,-1,-1,-1};
int ym[8]={-1,-1,0,1,1,1,0,-1};

#define xDim 128
#define yDim 128
#define popSize (xDim*yDim)

class Agent{
public:
	Agent *ancestor;
	int nrPointingAtMe;
	int born;
	double genome[GENES];
	int type;
	double fitness;
	Agent();
	~Agent();
	void setupRand(void);
	void setupSpecific(double Gx,double Gy);
	void inherit(Agent *from,double mutationRate);
	
    void LOD(FILE *F);
};

double PM[3][3]={{0.0,-1.0,1.0},{1.0,0.0,0.0},{1.0,0.0,0.0}};


vector<vector<Agent*>> population;

void play(Agent *A,Agent *B);
double euclidDist(double genomeA[GENES], double genomeB[GENES]);
double editDist(double genomeA[GENES], double genomeB[GENES]);

int main(int argc, const char * argv[]){
	int i,j,k,ni,nj;
	double maxFit;
    double minFit;
	char filename[1000];
    bool verbousOutput;
	FILE *pop;
	sprintf(filename,"POP_%s",argv[1]);
	mutationRate=atof(argv[2]);
	replacementRate=atof(argv[3]);
	double pA,pB,pC,pD;
	double Gx[3],Gy[3];
	pA=atof(argv[4]);
	pB=atof(argv[5]);
	pC=atof(argv[6]);
	pD=atof(argv[7]);
	for(int i=0;i<3;i++){
		Gx[i]=atof(argv[8+(i*2)]);
		Gy[i]=atof(argv[8+1+(i*2)]);
	}
    verbousOutput=(bool)atof(argv[14]);
    brh=(bool)atof(argv[15]);
	PM[0][1]=pA;
	PM[0][2]=pB;
	PM[1][0]=pC;
	PM[2][0]=pD;
	pop=fopen(filename,"w+t");
	srand((int)getpid()); // this is not cross-platofrm. Changed for condor compat.
	population.clear();
	population.resize(xDim);
	for(i=0;i<xDim;i++){
		population[i].resize(yDim);
		for(j=0;j<yDim;j++){
			Agent *A=new Agent;
			A->type=i%3;
			A->setupSpecific(Gx[A->type], Gy[A->type]);
			population[i][j]=A;
		}
	}
	int types[3]={0,0,0};
	for(globalUpdate=1;globalUpdate<generations;globalUpdate++){
		types[0]=types[1]=types[2]=0;
		//showPayoffs();
		maxFit=population[0][0]->fitness;
        minFit=population[0][0]->fitness;
        for(i=0;i<xDim;i++)
			for(j=0;j<yDim;j++){
				for(k=0;k<4;k++){
					play(population[i][j],population[(i+xm[2+k])&(xDim-1)][(j+ym[2+k])&(yDim-1)]);
				}
				types[population[i][j]->type]++;
			}
		for(i=0;i<xDim;i++)
			for(j=0;j<yDim;j++)
            {
				if(population[i][j]->fitness>maxFit)
					maxFit=population[i][j]->fitness;
                if(population[i][j]->fitness<minFit)
                    minFit=population[i][j]->fitness;
            }
        
		for(k=0;k<popSize*replacementRate;k++){
            //proportional selection
            if(maxFit-minFit>0.0){
                int counter = 0;
				do{
					i=rand()&(xDim-1);
					j=rand()&(yDim-1);
                    counter++;
				}while((randDouble>(population[i][j]->fitness-minFit)/(maxFit-minFit))&&(counter<xDim*yDim));
            }
			else {
				i=rand()&(xDim-1);
				j=rand()&(yDim-1);
			}
			do{
				ni=rand()&(xDim-1);
				nj=rand()&(yDim-1);
			}while((ni==i)&&(nj==j));
			population[ni][nj]->nrPointingAtMe--;
			if(population[ni][nj]->nrPointingAtMe<=0)
				delete population[ni][nj];
			population[ni][nj]=new Agent();
			population[ni][nj]->inherit(population[i][j], mutationRate);
		}
		fprintf(pop,"%i,%f,%f,%f",globalUpdate,(double)types[0]/(double)popSize,
				(double)types[1]/(double)popSize,
				(double)types[2]/(double)popSize);
        if((globalUpdate&1023)==1023 || globalUpdate==generations ||globalUpdate==1){
			printf("%i %f %f %f %f\n",globalUpdate,(double)types[0]/(double)popSize,
				   (double)types[1]/(double)popSize,
				   (double)types[2]/(double)popSize,maxFit);
            double refGenome[3][GENES];
            int count[3]={0,0,0};
            double mean[3],var[3];
            //set all to 0
            for(i=0;i<3;i++){
                mean[i]=0.0;
                var[i]=0.0;
                for(j=0;j<GENES;j++)
                    refGenome[i][j]=0.0;
            }
            //computer average genotype and mean of all genes used
            for(i=0;i<xDim;i++){
                for(j=0;j<yDim;j++){
                    count[population[i][j]->type]++;
                    for(int g=0;g<GENES;g++){
                        refGenome[population[i][j]->type][g]+=population[i][j]->genome[g];
                        mean[population[i][j]->type]+=population[i][j]->genome[g];
                    }
                }
            }
            //normalization step for the above
            for(i=0;i<3;i++){
                mean[i]/=(double)(count[i]*GENES);
                for(j=0;j<GENES;j++){
                    refGenome[i][j]/=(double)count[i];
                }
            }
            
            //compute actual variance for all three
            for(i=0;i<xDim;i++){
                for(j=0;j<yDim;j++){
                    for(int g=0;g<GENES;g++){
                        var[population[i][j]->type]+=pow(population[i][j]->genome[g]-mean[population[i][j]->type],2.0);
                    }
                }
            }
            for(int g=0;g<3;g++)
                var[g]/=count[g]*GENES;
            
            //now save all of it
            fprintf(pop,",%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                    editDist(refGenome[0],refGenome[1]),
                    editDist(refGenome[0],refGenome[2]),
                    editDist(refGenome[1],refGenome[2]),
                    euclidDist(refGenome[0],refGenome[1]),
                    euclidDist(refGenome[0],refGenome[2]),
                    euclidDist(refGenome[1],refGenome[2]),
                    var[0],var[1],var[2]
                    );
            //this is in here to still output the regular gen files just in case
            if(verbousOutput){
                //I removed a loop here, because it wasn't ding anything ... due to reused variable i
                char plants[1000];
                char herbs[1000];
                char micros[1000];

                sprintf(plants,"GEN_%i_TYPE_0_%s",globalUpdate,argv[1]);
                sprintf(herbs,"GEN_%i_TYPE_1_%s",globalUpdate,argv[1]);
                sprintf(micros,"GEN_%i_TYPE_2_%s",globalUpdate,argv[1]);
                FILE *P=fopen(plants,"w+t");
                
                FILE *H=fopen(herbs,"w+t");
                FILE *M=fopen(micros,"w+t");
                
                for(i=0;i<xDim;i++){
                    for(j=0;j<yDim;j++){
                        switch(population[i][j]->type){
                            case 0:
                                fprintf(P,"%f",population[i][j]->genome[0]);
                                fprintf(P,",%f",population[i][j]->genome[1]);
                                fprintf(P,"\n");
                                break;
                            case 1:
                                fprintf(H,"%f",population[i][j]->genome[0]);
                                fprintf(H,",%f",population[i][j]->genome[1]);
                                fprintf(H,"\n");
                                break;
                            case 2:
                                fprintf(M,"%f",population[i][j]->genome[0]);
                                fprintf(M,",%f",population[i][j]->genome[1]);
                                fprintf(M,"\n");
                                break;
                        }
                    }
                }
                fclose(P);
                fclose(H);
                fclose(M);
                
            }
        } else {
            fprintf(pop,",NA,NA,NA,NA,NA,NA,NA,NA,NA\n");
        }
	}
	fclose(pop);
	
	return 0;
}


Agent::Agent(){
	ancestor=NULL;
	nrPointingAtMe=1;
	born=globalUpdate;
	fitness=0.0;
}

Agent::~Agent(){
	if(ancestor!=NULL){
		ancestor->nrPointingAtMe--;
		if(ancestor->nrPointingAtMe==0)
			delete ancestor;
	}
}

void Agent::setupRand(void){
	int i;
	for(i=0;i<GENES;i++)
		genome[i]=randDouble;
}

void Agent::setupSpecific(double Gx,double Gy){
	genome[0]=Gx;
	genome[1]=Gy;
}


void Agent::inherit(Agent *from,double mutationRate){
	int i;
	from->nrPointingAtMe++;
	ancestor=from;
	type=from->type;
    for(i=0;i<GENES;i++){
        if(randDouble<mutationRate){
			genome[i]=randDouble;
        }
		else
			genome[i]=from->genome[i];
    }
    
}

void Agent::LOD(FILE *F){
	if(ancestor!=NULL)
		ancestor->LOD(F);
	else{
		fprintf(F,"generation");
		for(int i=0;i<GENES;i++)
			fprintf(F,",G%i",i);
		fprintf(F,"\n");
	}
	fprintf(F,"%i",born);
	for(int i=0;i<GENES;i++)
		fprintf(F,",%f",genome[i]);
	fprintf(F,"\n");
}




// *** play agents against each other

double euclidDist(double genomeA[GENES], double genomeB[GENES]){
    double M=0.0;
    for(int i=0;i<GENES;i++)
        M+=(genomeA[i]-genomeB[i])*(genomeA[i]-genomeB[i]);
    return pow(M,0.5);
}

double editDist(double genomeA[GENES], double genomeB[GENES]){
    double M=0.0;
    for(int i=0;i<GENES;i++)
        M+=pow((genomeA[i]-genomeB[i])*(genomeA[i]-genomeB[i]),0.5);
    return M;
}
void play(Agent *A,Agent *B){
    double M=editDist(A->genome,B->genome);
    if (!brh || A->type==2 || B->type==2){
        A->fitness+=PM[A->type][B->type]*(2-M);
        B->fitness+=PM[B->type][A->type]*(2-M);
    } else {
        A->fitness+=PM[A->type][B->type]*M;
        B->fitness+=PM[B->type][A->type]*M;
    }
    

}
