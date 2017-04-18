#ifndef FILTRE_H
#define FILTRE_H

//les valeurs min, max, moyen chaque sous_bande 1D
struct val_sous_bande_1D{

	double min;
	double max;
	double moyen;
};

//les valeurs moyen, variance et numéro de valeurs chaque sous_bande 2D
struct val_sous_bande_2D{
	
	int n;
	double moyen;
	double variance;
};

/*filtres de haar*/

void analyse_haar(double* x,int n);

void synthese_haar(double* x,int n);

/*filtres biorthogonaux 9/7*/

void analyse_97(double* x,int n);

void synthese_97(double* x,int n);

/*filtres biorthogonaux 9/7 lifting*/

void analyse_97_lifting(double* x,int n);

void synthese_97_lifting(double* x,int n);

/*analyse et reconstruit multirésolutions*/

void arm(double* x,int n,int niveau);

void iarm(double *x,int n,int niveau);

//calcul de la valeur minimale,maximal,moyenne d'une sous-bande issue d'une AMR
val_sous_bande_1D* calcul_sous_bande_1D(double *x, int n, int niveau);

//décomposition 2D
void analyse2D_97(double *m, int p, int ptotal);

/*analyse et reconstruit multirésolutions 2D d'une matrice carrée pxp, j : niveau*/

void amr2D_97(double* m,int p,int j);

void iamr2D_97(double* m,int p,int j);

/*ajuster après la amr2D_97 pour voir l'image*/

void ajuster_arm2D_97(double *m, int p,int j);

void iajuster_arm2D_97(double *m, int p,int j);

//calcul de la valeur moyenne et variance d'une sous-bande issue d'une AMR 2D
val_sous_bande_2D* calcul_sous_bande_2D(double* m, int p, int niveau);

//calcul des debits
double* calcul_debit(val_sous_bande_2D* vals, int niveau, int p, double b);

//quantifier image
//choice = true : utiliser fonction quantlm(x,n,nq)
//choice = false : utiliser fonction quantlm_idx(x,n,nq)
void quantifier(double* m, int p, double* debits, int niveau, bool choice);

//encoder sous une fichier binaire 'outs/lena_encode.bin'
void encoder(double* m, int p, double* debits, int niveau);

//calcul les erreurs par PSNR
double calcul_PSNR(double* m, double* n, int p);

void read_signal(double* x,int n,char* filename);

void save_signal(double* x,int n,char* filename);

#endif
