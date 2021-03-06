#include <stdio.h>
#include <stdlib.h>
#include "bmp.h"

double* charge_bmp256(char* fichier,int* largeur,int* hauteur) {
  FILE* fp;
  unsigned short bfType;
  unsigned long bfOffBits;
  unsigned long biWidth;
  unsigned long biHeight;
  unsigned short biBitCount;
  unsigned long biCompression;
  unsigned char* pixels;
  int pixelsSize;
  int x,y;
  double* m;

  fp=fopen(fichier,"rb");
  if (fp==NULL) {
    printf("charge_bmp256: impossible d'ouvrir le fichier %s en lecture !\n",fichier);
    return NULL;
  }

  // BMP specifications : http://members.fortunecity.com/shetzl/bmpffrmt.html
  // Lecture de l'ent�te
  fread(&bfType,sizeof(unsigned short),1,fp);
  if (bfType!=19778) {
    printf("charge_bmp256: le fichier %s n'est pas un fichier BMP !\n",fichier);
    fclose(fp);
    return NULL;
  }

  // Lecture de l'offset du debut du bitmap
  fseek(fp,10,SEEK_SET);
  fread(&bfOffBits,sizeof(unsigned long),1,fp);
  
  // Lecture de la largeur et de la hauteur
  fseek(fp,18,SEEK_SET);
  fread(&biWidth,sizeof(unsigned long),1,fp);
  *largeur=(int)biWidth;
  fread(&biHeight,sizeof(unsigned long),1,fp);
  *hauteur=(int)biHeight;

  // Verification que l'image est bien en mode 256 couleurs
  fseek(fp,28,SEEK_SET);
  fread(&biBitCount,sizeof(unsigned short),1,fp);
  if (biBitCount!=8) {
    printf("charge_bmp256: le fichier BMP %s n'est pas en mode 256 couleurs !\n",fichier);
    fclose(fp);
    return NULL;
  }

  // Verification que l'image n'est pas compress�e
  fseek(fp,30,SEEK_SET);
  fread(&biCompression,sizeof(unsigned long),1,fp);
  if (biCompression!=0) {
    printf("charge_bmp256: le fichier BMP %s est en mode compress� !\n",fichier);
    fclose(fp);
    return NULL;
  }

  // Allocation d'un bloc memoire pour lire les pixels et lecture de ceux-ci
  pixelsSize=(*largeur)*(*hauteur);
  pixels=(unsigned char *)malloc(pixelsSize);
  fseek(fp,bfOffBits,SEEK_SET);
  fread(pixels,pixelsSize,1,fp);

  // Copie dans un buffer de double et transposition des lignes
  m=(double *)calloc(pixelsSize,sizeof(double));
  for (y=0;y<*hauteur;y++) {
    for (x=0;x<*largeur;x++) {
      m[x+*largeur*(*hauteur-1-y)]=(double)pixels[x+*largeur*y];
    }
  }
  free(pixels);

  fclose(fp);

  return m;
}

int ecrit_bmp256(char* fichier,int largeur,int hauteur,double* m) {
  FILE* fp;
  unsigned short us;
  unsigned long ul;
  unsigned char uc;
  unsigned long i;
  int pixelsSize;
  int x,y;
  unsigned char* pixels;
  
  fp=fopen(fichier,"wb");
  if (fp==NULL) {
    printf("ecrit_bmp256: impossible d'ouvrir le fichier %s en �criture !\n",fichier);
    return 0;
  }

  pixelsSize=largeur*hauteur;

  // Conversion double => unsigned char
  pixels=(unsigned char *)malloc(pixelsSize);
  for (y=0;y<hauteur;y++) {
    for (x=0;x<largeur;x++) {
      double d;
      unsigned char c;
      d=m[x+largeur*y];
      if (d<0.0) c=0;
      else if (d>255.0) c=255;
      else c=(unsigned char)d;
      
      pixels[x+largeur*(hauteur-1-y)]=c;
    }
  }

  // Ecriture de l'ent�te standard
  // bfType
  us=19778;
  fwrite(&us,sizeof(unsigned short),1,fp);

  // bfSize 
  // taille image + taille BITMAPFILEHEADER + taille BITMAPINFOHEADER + taille palette
  ul=pixelsSize+14+40+256*4;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // bfReserved
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // bfOffBits
  // taille BITMAPFILEHEADER + taille BITMAPINFOHEADER + taille palette
  ul=14+40+256*4;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biSize
  ul=40;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biWidth
  ul=largeur;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biHeight
  ul=hauteur;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biPlanes
  us=1;
  fwrite(&us,sizeof(unsigned short),1,fp);

  // biBitCount
  us=8;
  fwrite(&us,sizeof(unsigned short),1,fp);

  // biCompression
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biSizeImage
  ul=pixelsSize;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biXPelsPerMeter
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biYPelsPerMeter
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biClrUsed
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // biClrImportant
  ul=0;
  fwrite(&ul,sizeof(unsigned long),1,fp);

  // Ecriture de la palette en niveaux de gris
  for (i=0;i<256;i++) {
    uc=i;
    fwrite(&uc,sizeof(unsigned char),1,fp);

    uc=i;
    fwrite(&uc,sizeof(unsigned char),1,fp);

    uc=i;
    fwrite(&uc,sizeof(unsigned char),1,fp);

    uc=0;
    fwrite(&uc,sizeof(unsigned char),1,fp);
  }

  // Ecriture de l'image
  fwrite(pixels,largeur*hauteur,1,fp);

  free(pixels);

  fclose(fp);
  return 1;
}

