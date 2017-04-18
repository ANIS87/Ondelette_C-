#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "filtre.h"
#include "bmp.h"


void prin(double x[],int n){

    for(int i=0;i<n;i++)

    printf("%2f\n ",x[i]);
    printf("\n");
}

void copy(double *x, double *y,int n){
	for(int i=0;i<n;i++) y[i] = x[i];
}


int main(){

	/* SIGNAL */	
	
	int n = 512;

	double* x=(double *)calloc(n,sizeof(double));

	read_signal(x,n,(char *)"test.txt");

	double* y = (double*)calloc(n,sizeof(double));
	copy(x,y,n);

	analyse_haar(x,n);
	synthese_haar(x,n);

	save_signal(x,n,(char *)"outs/resultat_haar.txt");

	copy(y,x,n);	
	analyse_97(x,n);
	synthese_97(x,n);

	save_signal(x,n,(char *)"outs/resultat_97.txt");

	copy(y,x,n);		
	analyse_97_lifting(x, n);
	synthese_97_lifting(x,n);

	save_signal(x,n,(char *)"outs/resultat_97_lift.txt");

	copy(y,x,n);		
	arm(x,n,2);
	//compter les valeurs min, max, moyen d'une sous_bande et sauver sur le fichier 'sous_bande.txt'
	val_sous_bande_1D* vals2 = calcul_sous_bande_1D(x, n, 2);	
	iarm(x,n,2);
	
	save_signal(x,n,(char *)"outs/resultat_arm_niveau2.txt");

	copy(y,x,n);	
	//9 = niveau maximal avec la longeur n = 512	
	arm(x,n,9); 
	val_sous_bande_1D* vals9 = calcul_sous_bande_1D(x, n, 9);
	iarm(x,n,9);
	
	save_signal(x,n,(char *)"outs/resultat_arm_niveau9.txt");

	
	/* IMAGE */

	char* fichier = (char *)"lena.bmp";
	char* fichier1 = (char *)"outs/lena1.bmp";
	char* fichier2 = (char *)"outs/lena2.bmp";

	int p;
	/*
	double* m = charge_bmp256(fichier,&p,&p);
	analyse2D_97(m,p));
	ecrit_bmp256(fichier1,p,p,m)
	*/
	//AMR 2D
	double* m = charge_bmp256(fichier,&p,&p);
	int niveau = 2;
	
	//decomposer image
	amr2D_97(m,p,niveau);

	//cacul les valeurs moyens et variances de chaque sous_bande
	val_sous_bande_2D* vals2d = calcul_sous_bande_2D(m, p, niveau);	

	//calcul les debits
	double* debits = calcul_debit(vals2d,niveau,p,2);

	//quantifier	
	quantifier(m, p, debits, niveau,true);

	//encoder	
	encoder(m,p,debits,niveau);

	//ajuster les valeurs pour voir l'image	
	ajuster_arm2D_97(m,p,niveau);
	

	ecrit_bmp256(fichier1,p,p,m);

	//reconstruire image
	iajuster_arm2D_97(m,p,niveau);
	iamr2D_97(m,p,niveau);

	ecrit_bmp256(fichier2,p,p,m);

	//calcul erreur PSNR
	double* origine = charge_bmp256(fichier,&p,&p);
	printf("\n PSNR = %f \n",calcul_PSNR(m, origine, p));

	return 0;
}


