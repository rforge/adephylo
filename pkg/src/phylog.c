#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Utils.h>
#include "adesub.h"

void gearymoran (int *param, double *data, double *bilis, 
    double *obs, double *result, double *obstot, double *restot);
    
void VarianceDecompInOrthoBasis (int *param, double *z, double *matvp,
    double *phylogram, double *phylo95,double *sig025, double *sig975,
    double *test1, double *test2, double*test3, double *test4, double *test5);

void MVarianceDecompInOrthoBasis (int *param, double *z, int *nvar, double *inertot, double *matvp,
    double *phylogram, double *phylo95,double *sig025, double *sig975,
    double *test1, double *test2, double*test3, double *test4, double *test5);

 
 void gearymoran (int *param, double *data, double *bilis, 
    double *obs, double *result, double *obstot, double *restot)
{
    /* Declarations des variables C locales */
    int nobs, nvar, nrepet, i, j, k, krepet, kvar ;
    int *numero;
    double provi;
    double *poili;
    double **mat, **tab, **tabperm;


    /* Allocation memoire pour les variables C locales */
    nobs = param[0];
    nvar = param [1];
    nrepet = param [2];
    vecalloc(&poili,nobs);
    taballoc(&mat,nobs,nobs);
    taballoc(&tab,nobs,nvar);
    taballoc(&tabperm,nobs,nvar);
    vecintalloc (&numero, nobs);

    /* D�finitions des variables C locales */
    k = 0;
    for (i=1; i<=nvar; i++) {
        for (j=1; j<=nobs; j++) {
            tab[j][i] = data[k] ;
            k = k+1 ;
       }
    }
    
    k = 0;
    provi = 0;
    for (j=1; j<=nobs; j++) {
        for (i=1; i<=nobs; i++) {
            mat[i][j] = bilis[k] ;
            provi = provi +  bilis[k];
            k = k+1 ;
       }
    }
    for (j=1; j<=nobs; j++) {
        for (i=1; i<=nobs; i++) {
            mat[i][j] = mat[i][j]/provi ;
       }
    }
    /* mat contient une distribution de fr�quence bivari�e */
    for (j=1; j<=nobs; j++) {
      poili[j] = 1/(double)(nobs);
    }
    /* poili contains uniform weights*/
    matmodifcn(tab,poili);
    /* data are standardized using uniform weights */
   
    for (kvar=1; kvar<=nvar; kvar++) {
        provi = 0;
        for (j=1; j<=nobs; j++) {
            for (i=1; i<=nobs; i++) {
                provi = provi + tab[i][kvar]*tab[j][kvar]*mat[i][j] ;
            }
        }
        obs[kvar-1] = provi;
    }
    k=0;
    /* les r�sultats se suivent par simulation */
    for (krepet=1; krepet<=nrepet; krepet++) {
        getpermutation (numero, krepet);
        matpermut (tab, numero, tabperm);
       
        for (kvar=1; kvar<=nvar; kvar++) {
            provi = 0;
            for (j=1; j<=nobs; j++) {
                for (i=1; i<=nobs; i++) {
                    provi = provi + tabperm[i][kvar]*tabperm[j][kvar]*mat[i][j] ;
                }
            }
            result[k] = provi;
            k = k+1;
        }
    }
    
    /* lib�ration m�moire locale */
    freevec(poili);
    freetab(mat);
    freeintvec(numero);
    freetab(tab);
    freetab(tabperm);
}
  

 void VarianceDecompInOrthoBasis (int *param, double *z, double *matvp,
    double *phylogram, double *phylo95,double *sig025, double *sig975,
    double *R2Max, double *SkR2k, double*Dmax, double *SCE, double *ratio)
{
    
    /* param contient 4 entiers : nobs le nombre de points, npro le nombre de vecteurs
    nrepet le nombre de permutations, posinega la nombre de vecteurs de la classe posi
    qui est nul si cette notion n'existe pas. Exemple : la base Bscores d'une phylog�nie a posinega = 0
    mais la base Ascores a posinega � prendre dans Adim
    z est un vecteur � nobs composantes de norme 1
    pour la pond�ration uniforme. matvp est une matrice nobsxnpro contenant en 
    colonnes des vecteurs orthonorm�s pour la pond�ration uniforme. En g�n�
    La proc�dure placera 
        dans phylogram les R2 de la d�composition de z dans la base matvp
        dans phylo95 les quantiles 0.95 des R2
        dans sig025 les quantiles 0.025 des R2 cumul�s
        dans sig975 les quantiles 0.975 des R2 cumul�s 
        
    Ecrit � l'origine pour les phylog�nies
    peut servir pour une base de vecteurs propres de voisinage */
        
    
    /* Declarations des variables C locales */
    int nobs, npro, nrepet, i, j, k, n1, n2, n3, n4;
    int irepet, posinega, *numero, *vecrepet;
    double **vecpro, *zperm, *znorm;
    double *locphylogram, *modelnul;
    double a1, provi, **simul, *copivec, *copicol;
    
   /* Allocation memoire pour les variables C locales */
    nobs = param[0];
    npro = param [1];
    nrepet = param [2];
    posinega = param[3];
    vecalloc (&znorm, nobs);
    vecalloc (&zperm, nobs);
    vecalloc (&copivec, npro);
    vecalloc (&copicol, nrepet);
    taballoc (&vecpro, nobs, npro);
    taballoc (&simul, nrepet, npro);
    vecalloc (&locphylogram, npro);
    vecalloc (&modelnul, npro);
    vecintalloc (&numero, nobs);
    vecintalloc (&vecrepet, nrepet);
    
    /* D�finitions des variables C locales */
    for (i = 1 ; i<= nobs; i++) znorm[i] = z[i-1];
    for (i = 1 ; i<= npro; i++) modelnul[i] = (double) i/ (double) npro;
    k = 0;
    for (j=1; j<=npro; j++) {
        for (i=1; i<=nobs; i++) {
            vecpro[i][j] = matvp[k] ;
             k = k+1 ;
       }
    }
    
   /* calcul du phylogramme observ� */
    for (j = 1; j<= npro; j++) {
        provi = 0;
        for (i=1; i<=nobs; i++)  provi = provi + vecpro[i][j]*znorm[i];
        provi = provi*provi/nobs/nobs;
        locphylogram[j] = provi;
   }
    for (i =1 ; i<= npro ; i++) phylogram[i-1] = locphylogram[i];
    /* calcul des simulations     
    Chaque ligne de simul est un phylogramme apr�s permutation des donn�es */
    
    for (irepet=1; irepet<=nrepet; irepet++) {
        getpermutation (numero, irepet);
        vecpermut (znorm, numero, zperm);
        provi = 0;
        for (j = 1; j<= npro; j++) {
            provi = 0;
            for (i=1; i<=nobs; i++)  provi = provi + vecpro[i][j]*zperm[i];
            provi = provi*provi/nobs/nobs;
            simul[irepet][j] = provi;
        }
    }
    /* calcul du test sur le max du phylogramme */
    for (irepet=1; irepet<=nrepet; irepet++) {
         for (j=1; j<=npro; j++) copivec[j] = simul[irepet][j];
         R2Max[irepet] = maxvec(copivec);
         provi=0;
         for (j=1; j<=npro; j++) provi = provi + j*simul[irepet][j];
         SkR2k[irepet] =provi;
         if (posinega>0) {
            provi=0;
            for (j=1; j<posinega; j++) provi = provi + simul[irepet][j];
            ratio[irepet] = provi;
        }
            
    }
    R2Max[0] = maxvec(locphylogram);
    provi=0;
    for (j=1; j<=npro; j++) provi = provi + j*locphylogram[j];
    SkR2k[0] =provi;
    if (posinega>0) {
            provi=0;
            for (j=1; j<posinega; j++) provi = provi + locphylogram[j];
            ratio[0] = provi;
   }
   /* quantiles 95 du sup */
    n1 = (int) floor (nrepet*0.95);
    n2 = (int) ceil (nrepet*0.95);
    for (i =1; i<=npro; i++) {
        for (irepet = 1; irepet<= nrepet; irepet++) {
            copicol[irepet] = simul [irepet][i];
        }
        trirap (copicol, vecrepet);
            phylo95[i-1] = 0.5*(copicol[n1]+copicol[n2]);
   }
   
  
  for (irepet=1; irepet<=nrepet; irepet++) {
        provi = 0;
        for (j=1; j<=npro; j++) {
            provi = provi + simul[irepet][j];
            copivec[j] = provi;
        }
        for (j=1; j<=npro; j++) simul[irepet][j] = copivec[j];
    } 
    n1 = (int) floor (nrepet*0.025);
    n2 = (int) ceil (nrepet*0.025);
    n3 = (int) floor (nrepet*0.975);
    n4 = (int) ceil (nrepet*0.975);
    /* quantiles 2.5 du cumul */
    for (i =1; i<=npro; i++) {
        for (irepet = 1; irepet<= nrepet; irepet++) {
            copicol[irepet] = simul [irepet][i];
        }
        trirap (copicol, vecrepet);
        sig025[i-1] = 0.5*(copicol[n1]+copicol[n2]);
        sig975[i-1] = 0.5*(copicol[n3]+copicol[n4]);
   }
   
    provi = 0;
    for (j=1; j<=npro; j++) {
        a1 = modelnul[j];
        provi = provi + locphylogram[j];
        locphylogram[j] = provi-a1;
        for (irepet = 1; irepet<= nrepet; irepet++) {
            simul [irepet][j] = simul [irepet][j]-a1;
        }
    }
    /* simul contient maintenant les cumul�s simul�s en �carts */
    /* locphylogram contient maintenant les cumul�s observ�s en �cart*/
    /* Dmax */
    for (j=1; j<=npro; j++) {
        for (irepet=1; irepet<=nrepet; irepet++) {
            for (j=1; j<=npro; j++) copivec[j] = simul[irepet][j];
            Dmax[irepet] = maxvec(copivec);
            provi=0;
            for (j=1; j<=npro; j++) provi = provi + copivec[j]* copivec[j];
            SCE[irepet] =provi;
        }
    }
    Dmax[0] = maxvec (locphylogram);
    provi=0;
    for (j=1; j<=npro; j++) provi = provi +locphylogram[j]*locphylogram[j];
    SCE[0] =provi;

   
   
    
    
    /* retour */
    
    freevec (znorm);
    freevec (modelnul);
    freevec(copivec);
    freevec(copicol);
    freevec (zperm);
    freetab (vecpro);
    freetab (simul);
    freevec (locphylogram);
    freeintvec (numero);
    freeintvec (vecrepet);
 }
    
void MVarianceDecompInOrthoBasis (int *param, double *z, int *nvar, double *inertot, double *matvp,
    double *phylogram, double *phylo95,double *sig025, double *sig975,
    double *R2Max, double *SkR2k, double*Dmax, double *SCE, double *ratio)
{
    
    /* param contient 4 entiers : nobs le nombre de points, npro le nombre de vecteurs
    nrepet le nombre de permutations, posinega la nombre de vecteurs de la classe posi
    qui est nul si cette notion n'existe pas. Exemple : la base Bscores d'une phylog�nie a posinega = 0
    mais la base Ascores a posinega � prendre dans Adim
    z est un vecteur � nobs composantes de norme 1
    pour la pond�ration uniforme. matvp est une matrice nobsxnpro contenant en 
    colonnes des vecteurs orthonorm�s pour la pond�ration uniforme. En g�n�
    La proc�dure placera 
        dans phylogram les R2 de la d�composition de z dans la base matvp
        dans phylo95 les quantiles 0.95 des R2
        dans sig025 les quantiles 0.025 des R2 cumul�s
        dans sig975 les quantiles 0.975 des R2 cumul�s 
        
    Ecrit � l'origine pour les phylog�nies
    peut servir pour une base de vecteurs propres de voisinage */
        
    
    /* Declarations des variables C locales */
    int nobs, npro, nrepet, i, j, k, n1, n2, n3, n4;
    int irepet, posinega, *numero, *vecrepet;
    double **vecpro, *zperm, **tabz;
    double *locphylogram, *modelnul;
    double a1, provi, **simul, *copivec, *copicol;
    
   /* Allocation memoire pour les variables C locales */
    nobs = param[0];
    npro = param [1];
    nrepet = param [2];
    posinega = param[3];
    vecalloc (&copivec, npro);
    vecalloc (&copicol, nrepet);
    taballoc (&vecpro, nobs, npro);
    taballoc (&tabz, nobs, nvar[0]);
    taballoc (&simul, nrepet, npro);
    vecalloc (&locphylogram, npro);
    vecalloc (&modelnul, npro);
   
    vecintalloc (&vecrepet, nrepet);
    
    /* D�finitions des variables C locales */
    
    for (i = 1 ; i<= npro; i++) modelnul[i] = (double) i/ (double) npro;
    k = 0;
    for (j=1; j<=npro; j++) {
        for (i=1; i<=nobs; i++) {
            vecpro[i][j] = matvp[k] ;
             k = k+1 ;
       }
    }
    k = 0;
    for (j=1; j<=nvar[0]; j++) {
        for (i=1; i<=nobs; i++) {
            tabz[i][j] = z[k] ;
             k = k+1 ;
       }
    }
    
   /* calcul du phylogramme observ� */
    
    
    for(k=1;k<=nvar[0];k++){     
      for (j = 1; j<= npro; j++) {
	provi = 0;
	for (i=1; i<=nobs; i++)  provi = provi + vecpro[i][j]*tabz[i][k];
	provi = provi*provi/nobs/nobs;
	locphylogram[j] =  locphylogram[j] + provi/inertot[0];
      }
    }
    for (i =1 ; i<= npro ; i++) phylogram[i-1] = locphylogram[i];
    /* calcul des simulations     
    Chaque ligne de simul est un phylogramme apr�s permutation des donn�es */
    
    for (irepet=1; irepet<=nrepet; irepet++) {
      aleapermutmat(tabz);
      provi = 0;
      for(k=1;k<=nvar[0];k++){     
	for (j = 1; j<= npro; j++) {
	  provi = 0;
	  for (i=1; i<=nobs; i++)  provi = provi + vecpro[i][j]*tabz[i][k];
	  provi = provi*provi/nobs/nobs;
	  simul[irepet][j] = simul[irepet][j]+provi/inertot[0];
	  
	}
      }
    }
    /* calcul du test sur le max du phylogramme */
    for (irepet=1; irepet<=nrepet; irepet++) {
         for (j=1; j<=npro; j++) copivec[j] = simul[irepet][j];
         R2Max[irepet] = maxvec(copivec);
         provi=0;
         for (j=1; j<=npro; j++) provi = provi + j*simul[irepet][j];
         SkR2k[irepet] =provi;
         if (posinega>0) {
            provi=0;
            for (j=1; j<posinega; j++) provi = provi + simul[irepet][j];
            ratio[irepet] = provi;
        }
            
    }
    R2Max[0] = maxvec(locphylogram);
    provi=0;
    for (j=1; j<=npro; j++) provi = provi + j*locphylogram[j];
    SkR2k[0] =provi;
    if (posinega>0) {
            provi=0;
            for (j=1; j<posinega; j++) provi = provi + locphylogram[j];
            ratio[0] = provi;
   }
   /* quantiles 95 du sup */
    n1 = (int) floor (nrepet*0.95);
    n2 = (int) ceil (nrepet*0.95);
    for (i =1; i<=npro; i++) {
        for (irepet = 1; irepet<= nrepet; irepet++) {
            copicol[irepet] = simul [irepet][i];
        }
        trirap (copicol, vecrepet);
            phylo95[i-1] = 0.5*(copicol[n1]+copicol[n2]);
   }
   
  
  for (irepet=1; irepet<=nrepet; irepet++) {
        provi = 0;
        for (j=1; j<=npro; j++) {
            provi = provi + simul[irepet][j];
            copivec[j] = provi;
        }
        for (j=1; j<=npro; j++) simul[irepet][j] = copivec[j];
    } 
    n1 = (int) floor (nrepet*0.025);
    n2 = (int) ceil (nrepet*0.025);
    n3 = (int) floor (nrepet*0.975);
    n4 = (int) ceil (nrepet*0.975);
    /* quantiles 2.5 du cumul */
    for (i =1; i<=npro; i++) {
        for (irepet = 1; irepet<= nrepet; irepet++) {
            copicol[irepet] = simul [irepet][i];
        }
        trirap (copicol, vecrepet);
        sig025[i-1] = 0.5*(copicol[n1]+copicol[n2]);
        sig975[i-1] = 0.5*(copicol[n3]+copicol[n4]);
   }
   
    provi = 0;
    for (j=1; j<=npro; j++) {
        a1 = modelnul[j];
        provi = provi + locphylogram[j];
        locphylogram[j] = provi-a1;
        for (irepet = 1; irepet<= nrepet; irepet++) {
            simul [irepet][j] = simul [irepet][j]-a1;
        }
    }
    /* simul contient maintenant les cumul�s simul�s en �carts */
    /* locphylogram contient maintenant les cumul�s observ�s en �cart*/
    /* Dmax */
    for (j=1; j<=npro; j++) {
        for (irepet=1; irepet<=nrepet; irepet++) {
            for (j=1; j<=npro; j++) copivec[j] = simul[irepet][j];
            Dmax[irepet] = maxvec(copivec);
            provi=0;
            for (j=1; j<=npro; j++) provi = provi + copivec[j]* copivec[j];
            SCE[irepet] =provi;
        }
    }
    Dmax[0] = maxvec (locphylogram);
    provi=0;
    for (j=1; j<=npro; j++) provi = provi +locphylogram[j]*locphylogram[j];
    SCE[0] =provi;

   
   
    
    
    /* retour */
    
    freetab (tabz);
    freevec (modelnul);
    freevec(copivec);
    freevec(copicol);
    
    freetab (vecpro);
    freetab (simul);
    freevec (locphylogram);

    freeintvec (vecrepet);
 }
    
