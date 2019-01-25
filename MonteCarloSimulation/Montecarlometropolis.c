#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char* argv[])
{

	int i, j, ii, jj, kk, tt;
	int N=40; //grid points in one direction
	int montecarlosteps=150e3;
	int saveafterthismanymontecarlosteps=100e3;
	int savefrequencybymod=100;
	int storageSizeRequired=(montecarlosteps-saveafterthismanymontecarlosteps)/savefrequencybymod;
		
	double V1 = 1;
	double V2 = -V1*2.0;
	double mu = atof(argv[1]);
	double TempConst = atof(argv[2]);
	int FileSaveNumber =atoi(argv[3]); //saves the data using this file number in RESULTS
	int FolderNumber=atoi(argv[4]); //folder number to save to

	double EEold=0.0; //Energy of whole system
	double EEnew=0.0; //Energy of whole system
	double EEStorage [storageSizeRequired];
	double theta; 
	double thetaStorage [storageSizeRequired];
	double EE2;
	double Cv [storageSizeRequired];
	double avgSigma;
	double avgSigmaStorage [storageSizeRequired];
		
	int StorageIndex = 0;

	long int seed;
	FILE* urand = fopen ( "/dev/urandom","r");
	fread(&seed,sizeof(long int),1,urand);
	fclose(urand);
	srand48(seed);

	int SS[N][N]; //stores spins

	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			double rn_temp = 0.5-drand48(); //generate random number centered about 0
			SS[i][j]=(int) ((rn_temp > 0.0) - (rn_temp < 0.0)); //this takes the sign of the random number
		}
	}

	/*
	FILE *fpa = fopen("OriginalSpin.out", "w");
	fwrite(SS, sizeof(int),N*N, fpa);
	fclose (fpa);
	
	printf("Original Spin Matrix \n");
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			printf("%d ", SS[i][j]);
		}
		printf("\n");
	} 
	*/
	
	//goes over all the points to find the energy
	for(i = 0; i < N; i++) {
		for(j = 0; j < N; j++) {
			EEold+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
					+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
							SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
							-mu*SS[i][j]);
		}
	}


for (tt=0; tt<montecarlosteps; tt++)
{
	
	for (kk=0; kk<10; kk++) //Try 30 flips 
	{
		EEnew=EEold;
		ii=(int) (drand48()*N); //choose random index 
		jj=(int) (drand48()*N); //choose random index
		
		//Substract the original spin contribution to the energy
		i=ii; j=jj;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=ii; j=((jj+1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=ii; j=((jj-1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=jj;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=jj;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=((jj+1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=((jj-1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=((jj+1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=((jj-1)%N+N)%N;
		EEnew-=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);

		SS[ii][jj]=-SS[ii][jj];	//flip it for calculation of new energy
		
		//Add the energy contribution of the flipped spin
		i=ii; j=jj;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=ii; j=((jj+1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=ii; j=((jj-1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=jj;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=jj;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=((jj+1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii-1)%N+N)%N; j=((jj-1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=((jj+1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		i=((ii+1)%N+N)%N; j=((jj-1)%N+N)%N;
		EEnew+=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
				+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
						SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
						-mu*SS[i][j]);
		
		SS[ii][jj]=-SS[ii][jj];	//flip it back because we need to check if the flip is appropirate
		
		if (EEnew-EEold<0 || drand48()<exp((EEold-EEnew)/TempConst)) //flip it because energy drops or because x<P
		{
			SS[ii][jj]=-SS[ii][jj]; //now actually flip it
			EEold=EEnew; //energy is accepted and updated to reflect that
		}
	}
		
	
		if (tt>=saveafterthismanymontecarlosteps &&  tt%savefrequencybymod==0)
		{	
		double ee = 0.0;
		avgSigma=0.0; 
		EE2=0.0;
		for(i = 0; i < N; i++) {
			for(j = 0; j < N; j++) {
				avgSigma+=SS[i][j]/((double) N*N);
				ee=(0.5*V1*SS[i][j]*(SS[i][((j+1)%N+N)%N]+SS[i][((j-1)%N+N)%N]+SS[((i-1)%N+N)%N][j]+SS[((i+1)%N+N)%N][j])
						+0.5*V2*SS[i][j]*(SS[((i-1)%N+N)%N][((j+1)%N+N)%N]+SS[((i-1)%N+N)%N][((j-1)%N+N)%N]+
								SS[((i+1)%N+N)%N][((j+1)%N+N)%N]+SS[((i+1)%N+N)%N][((j-1)%N+N)%N])
								-mu*SS[i][j]);
				EE2+=ee*ee;
			}
		}
		theta=0.5+avgSigma/2.0;
		thetaStorage[StorageIndex]=theta;
		EEStorage[StorageIndex]=EEold/((double) N*N);
		Cv[StorageIndex]=(EE2-EEStorage[StorageIndex]*EEStorage[StorageIndex])/(TempConst*TempConst);
		StorageIndex+=1;
		
		printf("%i mu:%f T:%f Iteration #%i. Theta:%f Avg. Cv:%f \n", StorageIndex, mu, TempConst, tt, theta, Cv[StorageIndex-1]);
		
			
		}
	 

}

	printf("\n");
	


	char buffer1[32];
	snprintf(buffer1, sizeof(char) * 32, "./RESULTS%i/ThetaValues%i.out",FolderNumber,FileSaveNumber);
	FILE *fp1 = fopen(buffer1, "wb");
	fwrite(thetaStorage, sizeof(double), storageSizeRequired, fp1);
	fclose (fp1);

	char buffer2[32];
	snprintf(buffer2, sizeof(char) * 32, "./RESULTS%i/Cv%i.out",FolderNumber,FileSaveNumber);
	FILE *fp2 = fopen(buffer2, "wb");
	fwrite(Cv, sizeof(double), storageSizeRequired, fp2);
	fclose (fp2);
 


return 0;
}
