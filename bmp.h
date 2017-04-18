#ifndef BMP_H
#define BMP_H

double* charge_bmp256(char* fichier,int* largeur,int* hauteur);

int ecrit_bmp256(char* fichier,int largeur,int hauteur,double* m);

#endif
