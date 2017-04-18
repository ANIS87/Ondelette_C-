#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "filtre.h"
#include "quantlm.h"

double arm2D_a,arm2D_b;

void init(double *x, int n){
	for(int i=0;i<n;i++) x[i] = 0;
}

void interpolation(double *x, double *y, int n){

	for(int i=0; i<n; i++) y[i] = 0;
	for(int i=0; i<n/2; i++){
		y[2*i] = x[i];
        y[2*i+1] = 0;
	}

}

void decimation(double *x, double *y, int n){

	for(int i=0; i<n; i++) y[i] = 0;
	for(int i=0; i<n/2; i++){
		y[i] = x[2*i];
	}

}

void filtrage(double  *x, int n, double *_h, int z, double *y){

    for(int i = 0; i < n; i++){
        for(int k = -z/2; k<=z/2;k++){

            if(i-k>=n) y[i] += _h[k+z/2]*x[2*(n-1) - (i-k)];
            else if(i-k<0) y[i] += _h[k+z/2]*x[k-i];
            else y[i] += _h[k+z/2]*x[i-k];

        }
    }
}

void analyse_haar(double* x,int n){

    double _h0[] = {1/sqrt(2), 1/sqrt(2), 0};
    double xb[n];
	init(xb,n);
    filtrage(x,n,_h0,3,xb);

    double xbd[n];
	init(xbd,n);
    decimation(xb,xbd,n);

    double _h1[] = {1/sqrt(2), -1/sqrt(2), 0};
    double xh[n];
	init(xh,n);
    filtrage(x,n,_h1,3,xh);

    double xhd[n];
	init(xhd,n);
    decimation(xh,xhd,n);

    for(int i=0;i<n/2;i++){
        x[i] = xbd[i];
        x[i+n/2] = xhd[i];
    }

}

void synthese_haar(double* x,int n){

    double xbd[n],xbh[n],xb[n],xg[n],xbi[n],xgi[n];
    double _g0[] = {0,1/sqrt(2),1/sqrt(2)};
    double _g1[] = {0,-1/sqrt(2),1/sqrt(2)};
	init(xbi,n);
    init(xgi,n);
    init(xbd,n);
    init(xbh,n);
    init(xb,n);
    init(xg,n);

    for(int i = 0;i<n/2;i++) {xbd[i] = x[i]; xbh[i] = x[i+n/2];}

    interpolation(xbd,xbi,n);

    filtrage(xbi,n,_g0,3,xb);

    interpolation(xbh,xgi,n);

    filtrage(xgi,n,_g1,3,xg);

    for(int i=0;i<n;i++) x[i] = xb[i] + xg[i];
}

void init_h0(double *_h0){
	_h0[0]=0.037828455507;
	_h0[1]=-0.023849465019;
	_h0[2]=-0.110624404418;
	_h0[3]=0.377402855613;
	_h0[4]=0.852698679009;
	_h0[5]=0.377402855613;
	_h0[6]=-0.110624404418;
	_h0[7]=-0.023849465019;
	_h0[8]=0.037828455507; 
}

void init_h1(double *_h1){
	_h1[0]=0.064538882629;
	_h1[1]=-0.040689417610;
	_h1[2]=-0.418092273222;
	_h1[3]=0.788485616406;
	_h1[4]=-0.418092273222;
	_h1[5]=-0.040689417610;
	_h1[6]=0.064538882629;
	_h1[7]=0.000000000000;
	_h1[8]=-0.000000000000;
}

void analyse_97(double* x,int n){
	
	double _h0[9];
	init_h0(_h0);

	double _h1[9];
	init_h1(_h1);

    double xb[n];
	init(xb,n);
    filtrage(x,n,_h0,9,xb);

    double xbd[n];
	init(xbd,n);
    decimation(xb,xbd,n);

    double xh[n];
	init(xh,n);
    filtrage(x,n,_h1,9,xh);

    double xhd[n];
	init(xhd,n);
    decimation(xh,xhd,n);

    for(int i=0;i<n/2;i++){
        x[i] = xbd[i];
        x[i+n/2] = xhd[i];
    }
}

void init_g0(double *_g0){
	_g0[0]=-0.064538882629;
	_g0[1]=-0.040689417610;
	_g0[2]=0.418092273222;
	_g0[3]=0.788485616406;
	_g0[4]=0.418092273222;
	_g0[5]=-0.040689417610;
	_g0[6]=-0.064538882629;
}

void init_g1(double *_g1){
	_g1[0]=0.000000000000;
	_g1[1]=-0.000000000000;
	_g1[2]=0.037828455507;
	_g1[3]=0.023849465019;
	_g1[4]=-0.110624404418;
	_g1[5]=-0.377402855613;
	_g1[6]=0.852698679009;
	_g1[7]=-0.377402855613;
	_g1[8]=-0.110624404418;
	_g1[9]=0.023849465019;
	_g1[10]=0.037828455507;
}

void synthese_97(double* x,int n){
	double xbd[n],xbh[n],xb[n],xg[n],xbi[n],xgi[n];
    double _g0[7];
	init_g0(_g0);
    double _g1[11];
	init_g1(_g1);
	init(xbi,n);
    init(xgi,n);
    init(xbd,n);
    init(xbh,n);
    init(xb,n);
    init(xg,n);

    for(int i = 0;i<n/2;i++) {xbd[i] = x[i]; xbh[i] = x[i+n/2];}

    interpolation(xbd,xbi,n);

    filtrage(xbi,n,_g0,7,xb);

    interpolation(xbh,xgi,n);

    filtrage(xgi,n,_g1,11,xg);

    for(int i=0;i<n;i++) x[i] = xb[i] + xg[i];
}

void analyse_97_lifting(double* x,int n){

	//x2N+1
	float a;
	a =-1.586134342;
	for(int i=1; i<n; i+=2) x[i] = a*x[i-1]+x[i]+a*x[i+1] ;
	x[n-1] = a*x[n-2]+x[n-1]+a*x[n-2] ;

	//X2N
	 a =-0.05298011854;
	for(int i=2; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1];
	x[0]=a*x[1]+x[0]+a*x[1];

	//X2N+1
	a =0.8829110762;
	for(int i=1; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1] ;
	x[n-1] = a*x[n-2]+x[n-1]+a*x[n-2] ;

	//X2N
	a=0.4435068522;
	for(int i=2; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1];
	x[0]=a*x[1]+x[0]+a*x[1];

	//echelle
	a=1/1.149604398;
	for (int i=0;i<n/2;i++){
		x[2*i] =x[2*i]/a ; 
		x[2*i+1] = a*x[2*i+1];
	}

	//copie 
	double z[n];
	for (int i=0;i<n/2;i++){
		z[i]=x[2*i];
		z[n/2+i]=x[2*i+1];
	}
	for (int i=0;i<n;i++)	x[i]=z[i];

}

void synthese_97_lifting(double* x,int n){

	//inverse 
	double z[n];
	for (int i=0;i<n/2;i++){
		z[2*i]=x[i];
		z[2*i+1]=x[n/2+i];
	}
	for (int i=0;i<n;i++){
	x[i]=z[i];

	} 

	//Mise à l'échelle
	float a;
	a=1.149604398;
	for (int i=0;i<n/2;i++){
		x[2*i] =x[2*i]/a ; 
		x[2*i+1] = a*x[2*i+1];
	}

	//Mise à jour"X2N
	a=-0.4435068522;
	for(int i=2; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1];
	x[0]=a*x[1]+x[0]+a*x[1];

	//"Prédiction" 2 X2N+1
	a=-0.8829110762;
	for(int i=1; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1] ;
	x[n-1] = a*x[n-2]+x[n-1]+a*x[n-2] ;

	//"Mise à jour" X2N
	a=0.05298011854;
	for(int i=2; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1];
	x[0]=a*x[1]+x[0]+a*x[1];

	//"Prédiction" 2 X2N+1
	a=1.586134342;
	for(int i=1; i<n; i+=2)	x[i] = a*x[i-1]+x[i]+a*x[i+1] ;
	x[n-1] = a*x[n-2]+x[n-1]+a*x[n-2] ;
}

void arm(double* x,int n,int niveau){
	
	for(int i=0;i<niveau;i++) analyse_97_lifting(x,n/(pow(2,i)));
	
}

void iarm(double *x,int n,int niveau){

	for(int i=niveau-1;i>=0;i--) synthese_97_lifting(x,n/(pow(2,i)));

}

void min_max_moyen(double *x, int n1, int n2, val_sous_bande_1D& val){
	val.min = val.max = x[n1];
	val.moyen = 0;
	for(int i=n1; i< n2; i++){
		val.moyen += x[i]/(n2-n1);		
		if(x[i]<val.min) val.min = x[i];
		if(x[i]>val.max) val.max = x[i];
	
	}
}

val_sous_bande_1D* calcul_sous_bande_1D(double *x, int n, int niveau){
	
	//val_sous_bande_1D vals[niveau+1];

	val_sous_bande_1D* vals=(val_sous_bande_1D *)calloc(niveau+1,sizeof(val_sous_bande_1D));

	int n1 = 0, n2;

	printf("\n les valeurs min max moyen des sous_bande 1D au niveau %d: \n",niveau);

	for(int i=niveau;i>=0;i--){
		
		n2 = n/pow(2,i);		
		min_max_moyen(x,n1,n2,vals[i]);
		printf("\n min max moyen sur (%d,%d) : %2f %2f %2f \n",n1,n2,vals[i].min,vals[i].max,vals[i].moyen);
		n1 = n2;

	}

	return vals;
}

void analyse2D_97(double *m, int p, int ptotal){

	for(int i = 0; i< p; i++){
		analyse_97_lifting(&m[i*ptotal],p);
	}
	
	double c[p];

	for(int j = 0; j< p; j++){
		for(int i = 0; i< p; i++) {c[i]= m[i*ptotal+j];}
		analyse_97_lifting(c,p);
		for(int i = 0; i< p; i++) {m[i*ptotal+j] = c[i];}
	}

}

void synthese2D_97(double *m, int p, int ptotal){
	
	double c[p];	
	
	for(int j = 0; j< p; j++){
		for(int i = 0; i< p; i++) {c[i]= m[i*ptotal+j];}
		synthese_97_lifting(c,p);
		for(int i = 0; i< p; i++) {m[i*ptotal+j] = c[i];}
	}	

	for(int i = 0; i< p; i++){
		synthese_97_lifting(&m[i*ptotal],p);
	}

}

void amr2D_97(double* m,int p,int j){
	
	for (int i = 0; i< j; i++) analyse2D_97(m, (int) (p/pow(2,i)), p);
	
}

void ajuster_arm2D_97(double *m, int p,int j){

	//trouver les min max du matrix d'approximation
	double min,max;
	min = max = m[0];
	for(int i=0;i<p/(pow(2,j));i++)
		for(int k=0;k<p/(pow(2,j));k++){
			if(m[i*p+k]<min) min = m[i*p+k];
			if(m[i*p+k]>max) max = m[i*p+k];
		}

	arm2D_a = 255/(max-min);
	arm2D_b = -1*arm2D_a*min;

	for(int i=0;i<p;i++)
		for(int k=0;k<p;k++){
			if( (i < p/pow(2,j)) && (k < p/pow(2,j)) ) 	m[i*p+k] = arm2D_a*m[i*p+k] + arm2D_b;
			else										m[i*p+k]+=127;
		}	
}

void iajuster_arm2D_97(double *m, int p,int j){
	
	for(int i=0;i<p;i++)
		for(int k=0;k<p;k++){
			if( (i < p/pow(2,j)) && (k < p/pow(2,j)) ) 	m[i*p+k] =(m[i*p+k] - arm2D_b)/arm2D_a;
			else										m[i*p+k]-=127;
		}	
}

void iamr2D_97(double* m,int p,int j){
	
	for (int i = j-1; i>=0; i--) synthese2D_97(m, (int) (p/pow(2,i)), p);

}

void moyen_variance(double* m, int l1, int l2, int c1, int c2, int ptotal, val_sous_bande_2D& val){

	double val1=0, val2=0;

	for(int i = l1; i < l2; i++)
		for(int j = c1; j < c2; j++){
			val1 += m[i*ptotal+j] * m[i*ptotal+j];
			val2 += m[i*ptotal+j];
		}

	val.n = (l2-l1)*(c2-c1);
	val1 /= val.n;
	val2 /= val.n;
	val.moyen = val2;
	val.variance = val1 - val2*val2;
	
}

val_sous_bande_2D* calcul_sous_bande_2D(double* m, int p, int niveau){

	int d=0,l1,l2,c1,c2,temp;
	//val_sous_bande_2D vals[3*niveau + 1];
	
	val_sous_bande_2D* vals=(val_sous_bande_2D *)calloc(3*niveau+1,sizeof(val_sous_bande_2D));

	//coef.approximation
	l1 = 0 ; l2 = p/pow(2,niveau); c1 = 0; c2 = p/pow(2,niveau);
	moyen_variance(m, l1, l2, c1, c2, p, vals[d]);
	d++;
	
	for(int k = niveau -1; k>= 0; k--){
		
		//coef.H
		temp = c1; c1 = c2; c2 = p/pow(2,k);
		moyen_variance(m, l1, l2, c1, c2, p, vals[d]);
		d++;

		//coef.V
		c2 = c1; c1 = temp; l1 = l2; l2 = p/pow(2,k);
		moyen_variance(m, l1, l2, c1, c2, p, vals[d]);
		d++;

		//coef.D
		c1 = c2; c2 = p/pow(2,k);
		moyen_variance(m, l1, l2, c1, c2, p, vals[d]);
		d++;
	}
	
	printf("\n les valeurs moyen, variance et numero d'element des sous_bande 2D au niveau %d: \n",niveau);

	for(int k = 0; k< (niveau*3+1); k++) printf("\n %f %f %d\n", vals[k].moyen, vals[k].variance,vals[k].n);

	return vals;
}

double* calcul_debit(val_sous_bande_2D* vals, int niveau, int p, double b){
	
	
	int n = 3*niveau + 1;
	int N = p*p;
	double* debits = (double *)calloc(n,sizeof(double));
	
	double val1 = 1;
	for(int j = 0; j<n; j++)
		val1 *= pow((vals[j].variance),(1.0*vals[j].n/N));

	//printf("\n %f --- %f \n", vals[0].variance*vals[0].variance,1.0*vals[0].n/N);

	for(int i = 0;i<n;i++)
		debits[i] = b + 0.5 * (log(vals[i].variance/val1)/log(2));	

	printf("\n les debits des sous_bande 2D au niveau %d: \n",niveau);

	for(int k = 0; k<n; k++) printf("\n %f \n", debits[k]);

	return debits;

}


//quantifier chaque sous bande
void quantifier_sous_bande(double* m, int l1, int l2, int c1, int c2, int p, double debit, bool choice){

	int n = (l2-l1)*(c2-c1);	
	double* x = (double*)calloc(n,sizeof(double));

	//copy valeur
	int d = 0;
	for(int i = l1; i<l2;i++)
	for(int j = c1; j<c2;j++)
		x[d++] = m[i*p+j];

	//calcul numero de quantifier
	int nq = floor(pow(2,debit));

	//quantifier	
	if(choice==true)
		quantlm(x,n,nq);
	else 
		quantlm_idx(x,n,nq);	
	//recopy
	d = 0;
	for(int i = l1; i<l2;i++)
	for(int j = c1; j<c2;j++)
		m[i*p+j] = x[d++];	
	
}

//quantifier une image
void quantifier(double* m, int p, double* debits, int niveau, bool choice){

	int d=0,l1,l2,c1,c2,temp;
	
	//coef.approximation
	l1 = 0 ; l2 = p/pow(2,niveau); c1 = 0; c2 = p/pow(2,niveau);
	quantifier_sous_bande(m, l1, l2, c1, c2, p, debits[d],choice);
	d++;

	for(int k = niveau -1; k>= 0; k--){
		
		//coef.H
		temp = c1; c1 = c2; c2 = p/pow(2,k);
		quantifier_sous_bande(m, l1, l2, c1, c2, p, debits[d],choice);
		d++;

		//coef.V
		c2 = c1; c1 = temp; l1 = l2; l2 = p/pow(2,k);
		quantifier_sous_bande(m, l1, l2, c1, c2, p, debits[d],choice);
		d++;

		//coef.D
		c1 = c2; c2 = p/pow(2,k);
		quantifier_sous_bande(m, l1, l2, c1, c2, p, debits[d],choice);
		d++;
	}
	
}

void encoder(double* m, int p, double* debits, int niveau){

	//encoder sous le fichier binaire
	FILE *fp=fopen("outs/lena_encode.bin","w");

	for (int i=0; i<p; i++)
	for (int j=0; j<p; j++) 
		fwrite(&m[i*p+j],sizeof(unsigned short),1,fp);
	
	fclose(fp);
	
}

double calcul_PSNR(double* m, double* n, int p){
	
	double eqm = 0;	
	for(int i = 0; i<p;i++)
	for(int j = 0; j<p;j++)
		eqm += (m[i*p+j] - n[i*p+j]) * (m[i*p+j] - n[i*p+j]);
	
	eqm /= (p*p);

	double psnr = 10*(log(255*255/eqm)/log(10));

	return psnr;
}

void read_signal(double* x,int n,char* filename){
	
	int i;
	FILE *fp=fopen(filename,"rt");
	for (i=0;i<n;i++) {
		float t;
		fscanf(fp,"%f\n",&t);
		x[i]=(double)t;
	}
	fclose(fp);
}

void save_signal(double* x,int n,char* filename){
	
	int i;
	FILE *fp=fopen(filename,"wt");
	for (i=0;i<n;i++) fprintf(fp,"%f\n",x[i]);
	fclose(fp);
}








