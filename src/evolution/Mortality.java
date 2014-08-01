/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package evolution;

/**
 *
 * @author timbrochier
 */

public class Mortality {

static float Mmax = 0.6f; // Mortalite naturelle max par jour
static float Mmax_compet = 0.5f; //Mortalite de competition max par jour
static float param_ajsutement_predation = 0.7f; // plus il est grand plus la mortalite diminue
static double M_predation, M_competition, M_senescence, M_peche;
static int Biom_landings;
public static double calcul_mortality(float fish_length, double fish_weigth, int fish_stade, int fish_age, float Carcapa, double effectif, float densite_biomasse_locale_TOTAL, float densite_biomasse_de_ce_banc_par_km2, float bathy,  float dt_jour) {

    // Taille des oeufs :
     float L_init = 0.1f; //cm

// 1) COMPETITION
// densite_biomasse_locale= somme des biomasses de ce bancs et des autres dans la maille de grille, divise par la surface de la maille en km2
// Carcapa = somme du carbonne disponible par km2 sur les 20 premiers metres de la colonne d'eau

if  (densite_biomasse_locale_TOTAL > 1.2*Carcapa){
    // (// 1.2 : on autorise un dépassement de 20% de la carcapa, sinon on risque
    // de faire des calculs de competition sur des effectif de surfplus très faible,
    // avec meme le risque de division par zero si round(effectif_surplus)=0

    //System.out.println(" COMPETITION! biomasse TOTALE par km2 DANS LA ZONNE = " + densite_biomasse_locale_TOTAL/1000000 + " tonnes ; densite_biomasse_de_ce_banc_par_km2 = " + densite_biomasse_de_ce_banc_par_km2/1000000 + " tonnes ; Carcapa = " + Carcapa/1000000 +" tonnes ");
//    double capamax; // effectif max pour une cohorte d'indiv de ce poids (fish_weigth)
    double surplus_biomasse=densite_biomasse_locale_TOTAL-Carcapa; // biomasse en surplus (qui depasse la Carcapa)

    double surplus_effectif=Math.round(surplus_biomasse/fish_weigth); // effectif en surplus (qui depasse la capamax)

    M_competition = mort_competition_function(surplus_effectif);
    effectif = Math.max(0,effectif-surplus_effectif) + ((int) Math.round(surplus_effectif*(1-Mmax_compet*M_competition*dt_jour)));
}

// 2) MORTALITE PAR PREDATION (SELON LA TAILLE):
// Correctif bathy : pour représenter le fait qu'il y a plus de prédateurs en haute mer, et + gros,
// on considere un facteur correctif à la taille du poisson avant calcul de la mortalite par predation :
float facteur_correctif_mort_bathy = Math.max(0.1f,1 - (float) Kinesis.predation_bathy(bathy));
float fish_length_corr = facteur_correctif_mort_bathy*fish_length;
M_predation = mort_size_function(fish_length_corr, fish_stade, fish_age)*dt_jour;

double beta = 1e-13;
M_senescence = beta*fish_age*fish_age*fish_age;
//if (bathy<-2000){
//System.out.println("bathy = " + bathy + " --> facteur_correctif_mort_bathy = " + facteur_correctif_mort_bathy);
//}


// 3) Mortalite de peche
// sur tout individu > 10cm qui se trouve sur le plateau continental
M_peche = mort_peche_function(fish_length)*dt_jour; //, bathy_actuelle)*dt_jour;

Biom_landings= (int) Math.round(effectif*(1 - Mmax*M_peche)*fish_weigth); // en grammes

effectif = (int) Math.round( (float) effectif*(1 - Mmax*(M_predation + M_senescence) - M_peche));
//            System.out.println(" nouvel_effectif = " + nouvel_effectif);
// System.out.println(" bathy = " + bathy + " ; facteur_correctif_mort_bathy = " + facteur_correctif_mort_bathy);
        // On prend un L_init qui correspond a la metamorphose = 60 jours
return effectif;
}

static double mort_peche_function(float fish_length){
    // il faudrait faire dependant du stade (capturabilite)
    // + il faut enregistrer les tonnes de captunes par zone
    double M_pech;

    if (fish_length>Simulation.Min_fishing_size){
        M_pech = Simulation.F_annuel / (float) Simulation.nb_jours_par_ans;
    }else
    {
         M_pech = 0;
    }
    return M_pech;
    };

static double mort_size_function(float fish_length, int stage, int age){
double M_pred;

    //float epsilon = L_init/10;
//    double M = (L_init-epsilon)/fish_length; // mortalite indicative, a multiplier par Mmax pou l'avoir par jour
    //double M = Math.exp((L_init-epsilon)/fish_length) / Math.exp((L_init-epsilon)/L_init); // mortalite indicative, a multiplier par Mmax pou l'avoir par jour

    // Okunishi et al 2012
if (stage == 0){ // Eggs
    M_pred = 0.57;
      }
else if (stage == 1){ // Yolk sac larvae
    M_pred = 0.3;
     }
else { // a partir du stade de feeding larvae, mortalite dependante de la taille : 
M_pred = 0.189*Math.exp(-fish_length/2.468);
}
    return M_pred;
}

static double mort_competition_function(double surplus_effectif){
double MM;

// Si competition, on applique une deuxieme fois la mortalite mort_size_function
        //M_competition = mort_size_function(fish_length,1); //% pour ne pas mettre la mortalite a L = Linit, trop forte
        //M_competition = M_predation; // Je considere la même mortalite

    // MORTALITE QUADRATIQUE :
//double alpha = 1/(capamax+1);
//double MM = Math.min(1,alpha*effectif);

    // MORTALITE EN SERIE f(n) = 1 - 1/(2^n)
    int n = (int) Math.round(surplus_effectif/1000);

     MM = 1 - 1/(Math.pow(2, n)) ;

    return MM;
}
}
