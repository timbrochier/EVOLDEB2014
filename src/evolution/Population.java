package evolution;

import java.io.*;

import java.util.*;
import java.io.IOException;
import util.OutpoutManager_netcdf;

public class Population {
// ------------------------------------------------------------------------------------------------------

    FileWriter POSITION_Writer;
    public static ArrayList Pop;
    public static ArrayList Pop_rec;
    public static int D_max_spatiotemp, D_max_spatial, D_max_temporel;
    // domaine de la population (en coordonnées de grille) :
    static double x_min, y_min, x_max, y_max;
    // taille des mailles (en Km)
    static float mean_grid_size;
    //float[] Env_mere;
    public static int nb_jours_ecoule_entre_roms_records;
    public static int[][][] successSerial;
    public static int[][] ageSuccessSerial;
    public static int[] pondeurMois;
    public static int superindiv_mort_total, compteur_id_poissons;
    boolean flagTrace = false;
    static int[][][] D_strates_spatiotemp;
    static int[][][] D_strates_spatiotemp_finale;
    int[] D_strates_temporelle;
    //int[][] D_strates_spatiale;
    public static int i_max, i_min, j_max, j_min, X, Y, sponge_grid;
    int[] liste_d;
    static int sortie_NORD, sortie_SUD, sortie_EST, sortie_OUEST;
    static float SEX_RATIO = 0.5f;
    OutpoutManager_netcdf Sim_writer;
    public static int[][] Biom_tot_2D, Biom_ichthyop_2D, Biom_petit_2D, Biom_moy_2D, Biom_grd_2D, Biom_landings_2D;
    public static float[] N_recmax_SI;
    public static int[] N_rec_SI;
    int nrec_jour;
    float nrec_max_jour;

    int nbzones;

    // ------------------------------------------------------------------------------------------------------

    public Population(OutpoutManager_netcdf bi) {
        this.Sim_writer = bi;
        Pop = new ArrayList();
        Pop_rec = new ArrayList();
        pondeurMois = new int[12];
        compteur_id_poissons = 0;
    }
// ------------------------------------------------------------------------------------------------------

    public void stepSerial() {
        System.out.println("Advection des larves et Migration des juveniles et adultes....");

        // CORRECTION BUG ADVECTION (19 Juin 2009) pour interpolation :
        nb_jours_ecoule_entre_roms_records = -1;

        int jour_init = 0;
        int jour_final = Simulation.nb_jours_par_ans;

// -----------------------  PARCOUR DES JOURS DE L'ANNEE  ----------------------
        for (int jour = jour_init; jour < jour_final; jour++) {

            // CHAQUE JOUR, Remise a zero des densite 2D : les premiers ne voient personne!
            mise_a_zero_densites_2D();

            //System.out.println("nb_jours_ecoule_entre_roms_records = " +nb_jours_ecoule_entre_roms_records);
            nb_jours_ecoule_entre_roms_records++;
            if ((jour % 10) == 0) {
                System.out.println("GEN " + Simulation.compteur_year + " :  jour de l'annee: " + jour + " -------------- Pop.size = " + Pop.size());
            }
            //***** DEBUT LOAD ENVIRONNEMENTAL DATA ----------------------------
            if (jour % Dataset_EVOL.dt_HyMo == 0) {
                try {
                    //                  if (Simulation.TEST_KINESIS && jour == 0) { //(jour % 30) == 0) { //
                    //                       CHAMPS FIXE PAR MOIS
//                        Dataset_EVOL.load_data(30*9, 1980); // <-- CHAMPS FIXE
                    //                    Dataset_EVOL.load_data(jour, Simulation.year);
                    //                  System.out.println("Lecture load_data load_data(1, 1980) ");
                    //            } else if (Simulation.TEST_KINESIS == false) {
                    //   System.out.println("Lecture load_data load_data(" + jour +", " + Simulation.year + ") ");
                    Dataset_EVOL.load_data(jour, Simulation.year);
//                    }
                } catch (IOException e) {
                    System.out.println("probleme load_data : " + e.getMessage());
                }
                //***** FIN LOAD ENVIRONNEMENTAL DATA ------------------------------
                // pour interpolation :
                nb_jours_ecoule_entre_roms_records = 0;
            }

            //**** ICI ON PARCOUR LA POPULATION --------------------------------
            for (int s = 0; s < Pop.size(); s++) {
//                    System.out.println("Poisson " + ss + " début...");
                Poissons Po = (Poissons) Pop.get(s);

                // Mise a zero des compteurs de reproduction : 
                Po.Nb_eggs_spawned_dt = 0;
                Po.Nb_SI_spawned_dt = 0;

                if ((jour >= Po.day_init || Po.age > 0)
                        && Po.living) {
//            System.out.println("Advection du poissons...." + ss );
                    Po.stepSerial(jour); // <-- ADVECTION OU DEPLACEMENT

                    // TEST DE RECRUTEMENT DES LARVES : REMPLACEE PAR LA SURVIE QUI DEPEND DE LA BOUFFE?
                    // * (a travers la mortalite size-dependant)
                    if (Po.age == Simulation.ageMinAtRecruitment) {

                        Po.isRetained_on_shelf = true;

                        if ((Po.lat> Dataset_EVOL.lat_max) || (Po.lat< Dataset_EVOL.lat_min)){
                            Po.living = false;
                            break;
                        }
                        //    Test_retention a 30 jours:
                        if (Simulation.flagRET_shelf) {
                            Po.isRetained_on_shelf = Math.abs(Po.bathyActuelle) < Dataset_EVOL.prof_talu;
                        }
                        int zone = (int) (Math.floor(Po.lat - Dataset_EVOL.lat_min));
// TEST : PROBA SUCCES DE REC (PAR COMETITION)

// ----------------------// COMPETITION AU RECRUTEMENT //-----------------------
                        
                        // Cette routine doit permettre de fixer le nombre max
                        // de Super Individu qui recrutent chaque jour
                        // ex : on permet au max un super-indiv qui recrute par 1/12eme de degré de latitude


//                        if ((Po.isRetained_on_shelf) && (N_rec_SI[zone] > N_recmax_SI[zone])) {
  //                          System.out.println(" COMPETITION AU RECRUTEMENT !!! - Before : zone =" + zone + "  N_recmax_SI[zone] =  " +  N_recmax_SI[zone] + " - N_rec_SI[zone] = " + N_rec_SI[zone] + "Po.living = "+Po.living);
    //                    }
                        int E_N_recmax_SI = (int) Math.floor(N_recmax_SI[zone]);
                        float reste = N_recmax_SI[zone] - E_N_recmax_SI;
                        int pow = 0;
                        if (Math.random() < reste) {
                        pow = 1;
                        }
                        int N_recmax_SI = (int) E_N_recmax_SI + pow;

                        if (Po.isRetained_on_shelf
                                && (N_rec_SI[zone] <= N_recmax_SI)
                                && (nrec_jour < nrec_max_jour))
                        {
                            N_rec_SI[zone] += 1;
                            nrec_jour +=1;
                            Po.living = true;
                            Po.isRecruited = true;
                       
                            try
                            {
                                Sim_writer.write_DataIndiv_NetCDF(Po);
                            }
                            catch(Exception ex)
                            {
                                System.out.println((new StringBuilder()).append("probleme dans ecriture netcdf DataIndiv : ").append(ex.getMessage()).toString());
                            }
                        } else {
                            Po.living = false;
                        }
                    }

                    if (Simulation.TEST_KINESIS == false) {
                        // PONTE ||
                        if (Po.DEB.adult_ready_to_spawn) {

                            if (Simulation.Strategie.equals("HE") && (Po.strates == null)) {
                                //    Repérage des strates
                                trouver_strat_ponte_potentielles(Po);
                            }
                            if (Simulation.Strategie.equals("HE") && (Po.first_spawn == false)) {
                                ponte_HomingEnv(Po, jour);
                            } else {
                                //System.out.println(" PONTE du poisson " + s );
                                ponte_Opp(Po, jour);
                            }
                        }
                    }

                    // ** POISSON MORT : SUPPRESSION DE LA LISTE (Sauvegarde? à voir +tard ...)
//                Po.living = Po.DEB. A VOIR
                    if (Po.living == false) {
//                    System.out.println("Poisson " + ss + " Po.living == false ");
                        Pop.remove(Po);
                    }
                } // if ((jour >= Po.day_init || Po.age > 0) && Po.living) {


                Biom_tot_2D[(int) Po.x][(int) Po.y] += (int) Po.S * Po.weight;
                if (Po.Body_length < 5) {
                    Biom_ichthyop_2D[(int) Po.x][(int) Po.y] += (int) Po.S * Po.weight;
                } else if ((Po.Body_length > 5) && (Po.Body_length <= 15)) {
                    Biom_petit_2D[(int) Po.x][(int) Po.y] += (int) Po.S * Po.weight;
                } else if ((Po.Body_length > 15) && (Po.Body_length <= 25)) {
                    Biom_moy_2D[(int) Po.x][(int) Po.y] += (int) Po.S * Po.weight;
                } else if (Po.Body_length > 25) {
                    Biom_grd_2D[(int) Po.x][(int) Po.y] += (int) Po.S * Po.weight;
                }


            } // fin de la boucle sur les individus

//+ AJOUTER LE TRUC "DEUXIEME PASSAGE" <--- ???
            if ((jour % 30) == 0) {
                System.out.println("Enregistrement netcdf");

            if (Simulation.outpout_netcdf_flag) {
                //System.out.println(" REMPLISSAGE DU NetCDF TEST  " );
                try {
                    Sim_writer.writeResults_NetCDF(Population.Pop);
                    Sim_writer.write_Biom2D_NetCDF();
                } catch (Exception ex) {
                    System.out.println("probleme dans creation netcdf 1 : " + ex.getMessage());
                }
            }
                // compteur de pas de temps pour le netcdf : 
            Simulation.compteur_t_rec_intra_year++;
                }
        }// fin de la boucle sur les jours
    }

// ------------------------------------------------------------------------------------------------------
    void ZoneGeo2Grid() {

        System.out.println("Zone de ponte initiale : ");
        System.out.println(" lat_min = " + Dataset_EVOL.lat_min + " ; lon_min = " + Dataset_EVOL.lon_min);
        System.out.println(" lat_max = " + Dataset_EVOL.lat_max + " ; lon_max = " + Dataset_EVOL.lon_max);

        double[] pGrid;
        pGrid = Dataset_EVOL.geo2Grid(Dataset_EVOL.lon_min, Dataset_EVOL.lat_min);
        x_min = pGrid[0];
        y_min = pGrid[1];

        pGrid = Dataset_EVOL.geo2Grid(Dataset_EVOL.lon_max, Dataset_EVOL.lat_max);
        x_max = pGrid[0];
        y_max = pGrid[1];

        System.out.println("x_min = " + Math.round(x_min) + " ; x_max = " + Math.round(x_max));
        System.out.println("y_min = " + Math.round(y_min) + " ; y_max = " + Math.round(y_max));

        // Nombre de strates totale :
        // Nombre de strates potentielles (ou h<3000m)
        i_max =
                (int) Math.round(x_max);
        i_min =
                (int) Math.round(x_min);
        j_max =
                (int) Math.round(y_max);
        j_min =
                (int) Math.round(y_min);




        int Nstrates_spatiale_total = 0;


        int Nstrates_spatiale_potentielles = 0;


        for (int i = i_min; i
                < i_max; i++) {
            for (int j = j_min; j
                    < j_max; j++) {
                if (Dataset_EVOL.isInWater(i, j)) {
                    Nstrates_spatiale_total++;


                    if (Dataset_EVOL.hRho[j][i] < Dataset_EVOL.prof_talu) {
                        //System.out.println("i = " + i + " ; j = " + j);
                        Nstrates_spatiale_potentielles++;


                    }
                }
            }
        }
        int Nstrates_spatiotemp_potentielles = Nstrates_spatiale_potentielles * Simulation.nb_strat_temporelles;

        System.out.println("Nstrates_spatiale_total = " + Nstrates_spatiale_total + " ; Nstrates_spatiale_potentielles = " + Nstrates_spatiale_potentielles);
        System.out.println("Nstrates_spatiotemp_potentielles = " + Nstrates_spatiotemp_potentielles);

        // Densite max :

        D_max_spatial = (int) Math.round((((double) 100 * Simulation.Pop_max) / ((double) Nstrates_spatiale_potentielles)));

        // D_max_temporel tq la population ne puisse pas être répartie sur moins du tiers de l'année (4 mois)
        D_max_temporel = (int) ((Simulation.Pop_max * Simulation.resolution_temporelle_scanEnv) / (30.f * 4.f));

        D_max_spatiotemp = (int) Math.round((((double) 100 * Simulation.Pop_max) / ((double) Nstrates_spatiotemp_potentielles)));
        System.out.println("D_max_spatiotemp : " + D_max_spatiotemp);
        System.out.println("D_max_spatial : " + D_max_spatial);
        System.out.println("D_max_temporel : " + D_max_temporel);
        System.out.println("flag_Densite_max : " + Simulation.flag_Densite_max);


        if ((D_max_spatiotemp < 1) && Simulation.flag_Densite_max) {
            System.out.println(" ATTENTION : Il faut fixer une Pop_max plus grande");
            System.exit(0);


        } // Taille des mailles
        float taille_zonale_domaine_lat_min = (float) util.Spheric_dist.distance_haversine(Dataset_EVOL.lat_min, Dataset_EVOL.lon_min, Dataset_EVOL.lat_min, Dataset_EVOL.lon_max);


        float taille_zonale_domaine_lat_max = (float) util.Spheric_dist.distance_haversine(Dataset_EVOL.lat_max, Dataset_EVOL.lon_min, Dataset_EVOL.lat_max, Dataset_EVOL.lon_max);


        float taille_latitudinale_domaine = (float) util.Spheric_dist.distance_haversine(Dataset_EVOL.lat_max, Dataset_EVOL.lon_min, Dataset_EVOL.lat_min, Dataset_EVOL.lon_min);
        //double taille_diagonale_domaine = util.Spheric_dist.distance_haversine(Dataset_EVOL.lat_min, Dataset_EVOL.lon_min, Dataset_EVOL.lat_max, Dataset_EVOL.lon_max);

        System.out.println("Taille du domaine : ");
        System.out.println("taille_zonale_domaine_lat_min : " + Math.round(taille_zonale_domaine_lat_min) + " Km");
        System.out.println("taille_zonale_domaine_lat_max : " + Math.round(taille_zonale_domaine_lat_max) + " Km");
        System.out.println("taille_latitudinale_domaine : " + Math.round(taille_latitudinale_domaine) + " Km");
        //System.out.println("taille_diagonale_domaine : " + taille_diagonale_domaine);
        System.out.println(" ");

        // taille des maille à lat_min :


        float taille_zonale_des_mailles_lat_min = taille_zonale_domaine_lat_min / (float) (x_max - x_min);


        float taille_zonale_des_mailles_lat_max = taille_zonale_domaine_lat_max / (float) (x_max - x_min);


        float taille_latitudinale_des_mailles = taille_latitudinale_domaine / (float) (y_max - y_min);

        System.out.println("taille_zonale_des_mailles_lat_min : " + taille_zonale_des_mailles_lat_min + " Km");
        System.out.println("taille_zonale_des_mailles_lat_max : " + taille_zonale_des_mailles_lat_max + " Km");
        System.out.println("taille_latitudinale_des_mailles : " + taille_latitudinale_des_mailles + " Km");

        mean_grid_size =
                (taille_zonale_des_mailles_lat_min
                + taille_zonale_des_mailles_lat_max
                + taille_latitudinale_des_mailles) / 3;
        System.out.println("taille_moyenne_des_mailles : " + mean_grid_size + " Km");

        // SPONGE : zone tampon autour de la zone de ponte, pour l'enregistrement des densite de position finale
        // Choisir le sponge tel que les individus au bord du domaine de ponte ne puisse pas sortir du sponge pendant leur advection
        // Conversion de sponge en nombre de mailles
        sponge_grid = (int) (Dataset_EVOL.sponge_km / mean_grid_size);
        System.out.println("Sponge = " + sponge_grid);

        // Taille du domaine utilié : X, Y
        X = (int) (1 + Population.x_max - Population.x_min);
        Y = (int) (1 + Population.y_max - Population.y_min);



        if (Simulation.flag_Densite_max || Simulation.flagGREG_dens) {
            //D_strates_spatiale = new int[X][Y];
            D_strates_temporelle = new int[Simulation.nb_strat_temporelles];
            D_strates_spatiotemp = new int[Simulation.nb_strat_temporelles][X][Y];


            if (Simulation.flagGREG_dens) {
                D_strates_spatiotemp_finale = new int[Simulation.nb_strat_temporelles][X + 2 * sponge_grid][Y + 2 * sponge_grid];


            }
            mise_a_zero_densites();
            System.out.println(" ");


        }
    }

//-------------------------
    void mise_a_zero_densites() {
        for (int d = 0; d
                < Simulation.nb_strat_temporelles; d++) {
            D_strates_temporelle[d] = 0;



            for (int i = 0; i
                    < X; i++) {
                for (int j = 0; j
                        < Y; j++) {
                    //D_strates_spatiale[i][j] = 0;
                    //for (int d = 0; d < Simulation.nb_strat_temporelles; d++) {
                    D_strates_spatiotemp[d][i][j] = 0;
                    //}


                }
            }

            if (Simulation.flagGREG_dens) {
                for (int i = 0; i
                        < (X + 2 * sponge_grid); i++) {
                    for (int j = 0; j
                            < (Y + 2 * sponge_grid); j++) {
                        D_strates_spatiotemp_finale[d][i][j] = 0;


                    }
                }
            }
        }
    }

    void mise_a_zero_densites_2D() {

nrec_jour = 0;
for(int zone=0; zone<nbzones; zone++){
    N_rec_SI[zone] = 0;
}

        for (int j = 0; j
                < Dataset_EVOL.ny; j++) {
            for (int i = 0; i
                    < Dataset_EVOL.nx; i++) {
                //D_strates_spatiale[i][j] = 0;
                //for (int d = 0; d < Simulation.nb_strat_temporelles; d++) {
                Biom_tot_2D[i][j] = 0;
                Biom_ichthyop_2D[i][j] = 0;
                Biom_petit_2D[i][j] = 0;
                Biom_moy_2D[i][j] = 0;
                Biom_grd_2D[i][j] = 0;
                Biom_landings_2D[i][j] = 0;
                //}


            }
        }
    }

    // Ponte initale tout les deux jours tout le long de l'année
    public void ponte_init() {
        int jour_debut_ponte = 0;


        int jour_fin_ponte = Simulation.nb_jours_par_ans;



        if (Simulation.TEST_KINESIS) {
            jour_debut_ponte = 0;
            jour_fin_ponte = jour_debut_ponte + 1;
        }
        
        Biom_tot_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];
        Biom_ichthyop_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];
        Biom_petit_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];
        Biom_moy_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];
        Biom_grd_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];
        Biom_landings_2D = new int[Dataset_EVOL.nx][Dataset_EVOL.ny];

        N_rec_SI = new int[nbzones];

        for (int jour = jour_debut_ponte; jour
                < jour_fin_ponte;
                jour++) {
            //System.out.println("Simulation.nb_jours_par_ans = " + Simulation.nb_jours_par_ans);
            //System.out.println("Simulation.nb_jours_par_ans/Simulation.resolution_temporelle_scanEnv =  " +Simulation.nb_jours_par_ans/Simulation.resolution_temporelle_scanEnv);
            for (int p = 0; p
                    < Simulation.nbFishesJour; p++) { // Ici "Simulation.nbFishes_pasdetemps" a été remplacé par "Simulation.nbFishesJour" le 21 mars 2013
                Poissons Po = new Poissons(jour);
//                compteur_id_poissons++;

                Pop.add(Po);
                //Biom_tot_2D[(int)Math.round(Po.x)][(int)Math.round(Po.y)]++;
            }
        }
    }

//------------------------------------------------------------------------------  
    public void copieListe() {

        // D'abord on supprime les poissons non recrutes de la liste Pop
        System.out.println("taille de la liste AVANT supression des non recrutes: " + Pop.size());



        int fin = Pop.size() - 1;


        for (int s = fin; s
                >= 0; s--) {
            Poissons Po = (Poissons) Pop.get(s);


            if (!Po.isRecruited) {
                Pop.remove(s);


            }

        }
        System.out.println("taille de la liste APRES supression des non recrutes: " + Pop.size());
        // Puis on passe (link) les elements de Pop a Pop_rec 
        Pop_rec = (ArrayList) Pop.clone();
        // enfin on efface Pop.	
        Pop.clear();

        // nombre d'indiv pondeur pour chaque mois : 
        // mise a zero :	



        for (int i = 0; i
                < 12; i++) {
            pondeurMois[i] = 0;


        }

        for (int p = 0; p
                < Pop_rec.size(); p++) {
            Poissons Po = (Poissons) Pop_rec.get(p);
            pondeurMois[Math.min(Math.round(Po.day_init / 30), 11)]++;


        }
    }

//------------------------------------------------------------------------------
// la ponte pour le Homming stricte,
// ponte par jour, les dates de pontes pouvant evoluer:
    public void ponte_HomingStrict() {
        System.out.println("PONTE HOMING STRICT... ");


        int pow, nb_oeuf, jourponte, day_init_temp, p, f, Nb_oeufs_pondu;


        boolean[] liste_indiv_ne_pouvant_pas_pondre;


        boolean tout_les_oeufs_sont_pondu, pondre;


        double patchiness_rad_maille, reste, nbponte_pour_cet_indiv; //, nb_nouveaux_super_indiv;
//        double patchiness_rad_maille = (patchiness_rad / Population.mean_grid_size);
        // nombre a pondre par chaque recrute
        // Jour de ponte: +- time_rad autour du jour de ponte de naissance
        // la nouvelle date Po.init sera ainsi la date de ponte pour l'année suivante...       

        // On parcour les individus survivants de l'année précédante (pop_rec) et on regarde leur date et lieux de naissance;
        // les nouveaux individus vont pondre dans un rayon autour de ce lieu et date
        Nb_oeufs_pondu = 0;
        liste_indiv_ne_pouvant_pas_pondre = new boolean[Pop_rec.size()];

        //tout_les_oeufs_sont_pondu = false;
        //int rec_dens = Po.isRecruitedGREGdens ? 1 : 0;



        for (p = 0; p
                < Pop_rec.size(); p++) {
            liste_indiv_ne_pouvant_pas_pondre[p] = false;


        }

        while (Nb_oeufs_pondu < Simulation.Pop_max) {

            int nb_oeufs_restants_a_pondre = Simulation.Pop_max - Nb_oeufs_pondu;


            int nb_poisson_pouvant_pondre = 0;


            for (p = 0; p
                    < Pop_rec.size(); p++) {
                if (!liste_indiv_ne_pouvant_pas_pondre[p]) {
                    nb_poisson_pouvant_pondre++;


                }
            }

            for (p = 0; p
                    < Pop_rec.size(); p++) {
                if (liste_indiv_ne_pouvant_pas_pondre[p] == false) {
                    //System.out.println("Poisson " + p);
                    Poissons Po = (Poissons) Pop_rec.get(p);
                    pondre = true;

                    patchiness_rad_maille =
                            ((double) Po.rayon_exploration_spatial / (double) Population.mean_grid_size);

                    // nombre d'oeufs que va pondre ce poisson:

                    // ! BUG QUI FAIT PONDRE PLUS EN DEBUT D'ANNEE :
                    //ERROR    //nbponte_pour_cet_indiv = (double) 2 * (((float) Simulation.Pop_max - (float) Nb_oeufs_pondu) / ((float) Pop_rec.size() - (float) p));
                    // 1er correction du BUG ( revoir pour le cas ou beaucoup de poissons rec ne pondent pas)
                    //nbponte_pour_cet_indiv = (double) ((float) Simulation.Pop_max / ((float) Pop_rec.size()));
                    // 2eme correction du BUG
                    nbponte_pour_cet_indiv = (double) ((float) nb_oeufs_restants_a_pondre / ((float) nb_poisson_pouvant_pondre));

                    reste = nbponte_pour_cet_indiv - Math.floor(nbponte_pour_cet_indiv);
                    pow = 0;


                    if (Math.random() < reste) {
                        pow = 1;


                    }
                    nb_oeuf = (int) (Math.floor(nbponte_pour_cet_indiv) + pow);

                    // Ponte effective des oeufs :

                    jourponte = get_jourponte(Po.day_init, Math.round(Math.round(Po.rayon_exploration_temporel)));

                    // CAS DE DENSITE_MAX :


                    if (Simulation.flag_Densite_max) {
                        int s = 0;


                        while (D_strates_temporelle[jourponte] > D_max_temporel) {
                            jourponte = get_jourponte(Po.day_init, Math.round(Math.round(Po.rayon_exploration_temporel)));
                            s++;

                            if (s > 7000) {
                                System.out.println("Rayon d'explo temporelle trop petit pour sortir de la periode de forte densite...  : ");


                                break;


                            }
                        }

                        int ii = Math.round(Math.round(Po.x_init));


                        int jj = Math.round(Math.round(Po.y_init));


                        int i = ii - Population.i_min;


                        int j = jj - Population.j_min;


                        int d = Math.round(jourponte / Simulation.resolution_temporelle_scanEnv);
                        pondre = (D_strates_spatiotemp[d][i][j] < D_max_spatiotemp)
                                && (D_strates_temporelle[jourponte] < D_max_temporel); // pondre si Densite < Densite max



                        if (pondre) {
                            // ON INCREMENTE LA MATRICE DES DENSITEES :
                            D_strates_spatiotemp[d][i][j] += nb_oeuf;
                            D_strates_temporelle[jourponte] += nb_oeuf;



                        } else {
                            tout_les_oeufs_sont_pondu = false;
                            liste_indiv_ne_pouvant_pas_pondre[p] = true;


                        }
                    }

                    if (pondre) {
                        // POUR VRAIMENT EMPECHER LA PONTE DANS LES ZONES OU LA DENSITE EST TROP FORTE IL FAUDRAIT FAIRE LE WHILE OUTZONE ICI.
                        for (f = 0; f
                                < nb_oeuf; f++) {
                            Poissons indiv = new Poissons(Po.x_init, Po.y_init, Po.depth_init, jourponte, patchiness_rad_maille, Po.Tolerances_spatiotemps, Po.Tolerances_TS);
//                            compteur_id_poissons++;

                            Pop.add(indiv);
                            Nb_oeufs_pondu++;

                        }


                    }
                }

            }
        }
        // FIN PONTE


        // Merci Jeremie pour la suite:
        // 1-on classe Po_rec dans l'ordre des nouvelles dates de naissance:
        Collections.sort(Pop, new MonComparator());
        // 2- on initialise les données environementale
        day_init_temp = -1;


        for (p = 0; p
                < Pop.size(); p++) {
            Poissons Po = (Poissons) Pop.get(p);


            if (Simulation.record_initial_env && (Po.day_init > day_init_temp
                    && (Po.day_init % Dataset_EVOL.dt_HyMo == 0))) {
                try {
                    Dataset_EVOL.load_data(Po.day_init, Simulation.year);


                } catch (IOException e) {
                    System.out.println("probleme dans ponte_init : " + e.getMessage());


                }
            }
            //System.out.println("init()...");
            Po.init(Po.day_init);
            day_init_temp = Po.day_init;


        }
    }

//------------------------------------------------------------------------------
    public static int get_jourponte(int birthday, int time_rad) {
        int decalage = Math.round(Math.round(Math.random() * time_rad * 2));


        int jourponte = birthday + decalage - time_rad;


        if (jourponte < 0) {
            jourponte = jourponte + Simulation.nb_jours_par_ans;


        }
        if (jourponte >= Simulation.nb_jours_par_ans) {
            jourponte = jourponte - Simulation.nb_jours_par_ans;


        }
        return (jourponte);







    }

// Pour le classement de la pop dans l'ordre de naissance:
    public class MonComparator implements Comparator {
// (de Jeremie)
//    MonComparator() {
//    }

        public int compare(Object arg0, Object arg1) {
            Poissons p1 = (Poissons) arg0;
            Poissons p2 = (Poissons) arg1;
            if (p1 == null) {
                return 1;
            }
            if (p2 == null) {
                return -1;
            }
            int result = p1.day_init - p2.day_init;
            return result;
        }
    }

// la ponte pour le Homing Environemental:
    public void ponte_HomingEnv(Poissons Po, int day) {    // (int d, int X, int Y, double patchiness) {

//        System.out.println("Début ponte HE....");

        int d;



        int i, j, k, f, pow, nb_oeuf;
//        int jourponte = d * Simulation.resolution_temporelle_scanEnv;


        double proba_ponte_strate, reste;
//       % double proba_ponte = ((double) (Simulation.nbFishesJour * Simulation.dt_jour) / ((double) Simulation.Nb_strates_date[d]));

        // La probabilité que le poisson p ponde en chaque position potentielle est :
        proba_ponte_strate = 10; //Simulation.nbponte_par_indiv / Po.Nb_strates;

//        System.out.println("nbponte_par_indiv  = " + Simulation.nbponte_par_indiv);
//        System.out.println("Po.Nb_strates = " + Po.Nb_strates);
//        System.out.println("N_ponte_pour_fermeture_pop = " + N_ponte_pour_fermeture_pop);

        reste = proba_ponte_strate - Math.floor(proba_ponte_strate);

// TECHNIQUE DU TABLEA 4D :

// au lieu de parcourir tout le tableau en entier comme dans EVOL, on regarde
        // directement si la position du poisson est ok pour ponte :

        d = Math.round(day / Simulation.resolution_temporelle_scanEnv);
        i = Math.round(Math.round(Po.x)) - i_min;
        j = Math.round(Math.round(Po.y)) - j_min;
        k = 0; // PAS DE HOMING SUR VERTICALE



        double fx, fy;
        fx = Math.max(Po.x, x_min);
        fy = Math.max(Po.y, y_min);

        fx = Math.min(fx, x_max);
        fy = Math.min(fy, y_max);

        /*
        for (d = 0; d
        < Simulation.nb_strat_temporelles; d++) {
        for (i = 0; i
        < X; i++) {
        for (j = 0; j
        < Y; j++) {
        for (k = 0; k
        < Dataset_EVOL.prof_potentielles.length; k++) {
         *
         */
// Si on est dans la zone de ponte potentielle (X, Y)


        if ((i >= 0) && (j >= 0) && (i < i_max - 1) && (j < j_max - 1)) {
            // Si pour ce poisson, cette date a cet endroit on a un envOK...:
            if ((Po.strates[d][i][j]) == true) {

                System.out.println("POOOOOONNNNTE HE !!!");

                // le jour de ponte est :


                int jourponte = d * Simulation.resolution_temporelle_scanEnv;

                pow = 0;


                if (Math.random() < reste) {
                    pow = 1;


                }
                nb_oeuf = (int) (Math.floor(proba_ponte_strate) + pow);



                if (Simulation.flag_Densite_max || Simulation.flagGREG_dens) {
                    D_strates_spatiotemp[d][i][j] += nb_oeuf;
                    D_strates_temporelle[d] += nb_oeuf;


                }

                for (f = 0; f
                        < nb_oeuf; f++) {
//                                jourponte = get_jourponte(jour, Simulation.alea_temps);
                    // rayon de patchiness = 0,5 mailles pour assurer une ponte sur l'ensemble du territoire
                    Poissons indiv = new Poissons(fx, fy, Dataset_EVOL.prof_potentielles[k], jourponte, 0.5, Po.Tolerances_spatiotemps, Po.Tolerances_TS);
//                    compteur_id_poissons++;

                    Pop.add(indiv);
                    Simulation.Nb_oeufs_pondu++;


                }
            } else {
                System.out.println("..envir pas bon.");


            }
        } else {
            //         System.out.println("..condition 1 pas bonne : out of domain");
        }
    }
//------------------------------------------------------------------------------
// Ponte pour Opportunisme : 

    public void ponte_Opp(Poissons Po, int jourponte) {

        // EVOL-DEB : Ponte_opp = ponte si Po.ready_to_spawn et Env_ok
        // (Env_ok = independament de Env de naissance)
        int i, j, k, f, pow, S_newSI;
        long nb_nouveaux_super_indiv, NEggs;


        double N_ponte_pour_fermeture_pop, reste;

        // 1 - INDIVIDU PRET A PONDRE? (double check)


        if (Po.DEB.adult_ready_to_spawn) {
            // System.out.println("prêt à pondre...");

            // 2 - ENVIRIONNEMENT OK POUR PONDRE?
            // Passage en grille "scannee" (pour lectures rapide des facteurs environnement)
            i = Math.round(Math.round(Po.x)) - i_min;
            j = Math.round(Math.round(Po.y)) - j_min;


            int d = Math.round(jourponte / Simulation.resolution_temporelle_scanEnv);
            k = 0; // PAS DE HOMING SUR VERTICALE

            // COMPAR ENVIRONNEMENT ICI AVEC CONDITIONS ENVIR POUR PONTE
            // (ON PEUT OPTIMISER EN NE FAISANT CES CALCULS QU'UNE FOIS (min - max)
            //double waterdensity_min = Simulation.DensOpp - Simulation.dens_marge;
            //double waterdensity_max = Simulation.DensOpp + Simulation.dens_marge;
            //double waterdensity = util.BuoyancyScheme.waterDensity(Po.salt, Po.temp);


            double temperature_min = Simulation.Topp - Simulation.temp_marge;


            double temperature_max = Simulation.Topp + Simulation.temp_marge;


            double salt_min = Simulation.Sopp - Simulation.salt_marge;


            double salt_max = Simulation.Sopp + Simulation.salt_marge;



            boolean env_ok;

            env_ok = ((Po.temp > temperature_min && Po.temp < temperature_max) && // 2 HOMING TEMPERATURE
                    (Po.salt > salt_min && Po.salt < salt_max));// && // 3 HOMING SALINITE
            //(Th_opp > Th[d][i][j] - Th_marge && Th_opp < Th[d][i][j] + Th_marge) && // 5 HOMING PROF THERMO
            //(Po.bathyActuelle > Simulation.bathymin_opp && Po.bathyActuelle < Simulation.bathymax_opp); // && // 4 HOMING BATHYMETRY
            //(waterdensity > waterdensity_min && waterdensity < waterdensity_max) &&



            if (env_ok) {

                // A) CREATION DES NOUVEAUX SUPER-INDIVIDUS --------------------

                double fx, fy;
                fx = Po.x;
                fy = Po.y;

                if ((fx < Population.x_min) || // ajout 30 septembre; 5 Nov: modif "Population.sponge_grid"
                        (fx > Population.x_max)
                        || (fy < Population.y_min) || // ajout 30 septembre
                        (fy > Population.y_max)) {
                    //System.out.println("PONTE ANNULEE CAR HORS DOMAINE : x = " + Math.round(fx) + " , y = " + Math.round(fy) + " lon = " + Po.lon + " ; lat = " + Po.lat);
                } else if (Dataset_EVOL.isCloseToCost(fx, fy)) {
                    //System.out.println("PONTE ANNULEE CAR CLOSE TO COAST : x = " + fx + " , y = " + fy);
                    //System.out.println("(Dataset_EVOL.isCloseToCost(fx, fy)) = " + (Dataset_EVOL.isCloseToCost(fx, fy)));
                } else {
                    //System.out.println("on est bon pour pondre ....  = " );

                    // B/ CALCUL DU NOMBRE D'OEUF ET DE SUPER-INDIVIDUS ------------
                    Po.DEB.spawn(); // calcule Nb_oeuf et met a zero le reproduction Buffer (Er)
                    // Le nombre de eggs pondu par le super indiv poisson Po est :
                    NEggs = Math.round(((float) Po.S) * ((float) Po.DEB.Nb_oeuf) * SEX_RATIO);

                    // Le nombre des nouveau super individu pondu est calculé de façon a conserver un nombre total de super-individus proche de Pop_max :
                /*N_ponte_pour_fermeture_pop = Math.max( (double) (Simulation.Pop_max / Pop.size()), 1);
                    reste = N_ponte_pour_fermeture_pop - Math.floor(N_ponte_pour_fermeture_pop);
                    pow = 0;
                    if (Math.random() < reste) {
                    pow = 1;
                    }
                    nb_nouveaux_super_indiv = (int) (Math.floor(N_ponte_pour_fermeture_pop) + pow);
                     * */

                    nb_nouveaux_super_indiv = (int) Math.max(Math.floor(NEggs / Simulation.nbeggs_max_per_SI), 1);

//        System.out.println("nbponte_par_indiv  = " + Simulation.nbponte_par_indiv);
//        System.out.println("Po.Nb_strates = " + Po.Nb_strates);
//        System.out.println("N_ponte_pour_fermeture_pop = " + N_ponte_pour_fermeture_pop);

                    // L'effectif S des nouveau super individu pondu est :
                    S_newSI = (int) Math.round(NEggs / nb_nouveaux_super_indiv);

                    //      System.out.println("Mise à zero du buffer (E_R) et ponte de " + NEggs + " oeufs, repartits en " + nb_nouveaux_super_indiv + " super indiv");

// debug


                    if (NEggs < 0) {
                        System.out.println("NEggs = " + NEggs + " ; Po.S = " + Po.S + " ; Po.DEB.Nb_oeuf = " + Po.DEB.Nb_oeuf);


                    }

                    for (f = 0; f
                            < nb_nouveaux_super_indiv; f++) {

                        // rayon de patchiness = 0,5 mailles pour assurer une ponte sur l'ensemble du territoire
                        Poissons indiv = new Poissons(fx, fy, Dataset_EVOL.prof_potentielles[k], jourponte, 1.f, S_newSI);
//                        compteur_id_poissons++;

                        Pop.add(indiv);
                        //  System.out.println("POOOOOONNNNTE OPORTUNISTE!!!");
                        //               Simulation.Nb_oeufs_pondu++;


                    }
                    Po.first_spawn = false;

                    Po.Nb_eggs_spawned_dt = NEggs;
                    Po.Nb_SI_spawned_dt = (int) nb_nouveaux_super_indiv;

                    //  System.out.println("NEggs = " + NEggs + " ; nb_nouveaux_super_indiv = " + nb_nouveaux_super_indiv );
                    //  System.out.println("Po.Nb_eggs_spawned_dt = " + Po.Nb_eggs_spawned_dt + " ; Po.Nb_SI_spawned_dt = " + Po.Nb_SI_spawned_dt );


                    // On enregistre la densite de ponte en ce point (A METTRE AU POINT)


                    if (Simulation.flag_Densite_max || Simulation.flagGREG_dens) {
                        D_strates_spatiotemp[d][i][j] += nb_nouveaux_super_indiv;
                        D_strates_temporelle[d] += nb_nouveaux_super_indiv;


                    }
                }
            } else {
//((Po.temp > temperature_min && Po.temp < temperature_max) && // 2 HOMING TEMPERATURE
                // System.out.println("... mais environnement pas acceptable = Po.temp = " + Po.temp  + " ; Po.salt = " + Po.salt);
            }
        }
    }

    // ----------------------------------------------------------------------------
    void pop_init_HE() {
        //int day_init_temp;

        // Pour ne pas enregistrée les données envir mais garder celle des parents..
        // Simulation.record_initial_env = false;

        // Merci Jeremie pour la suite:     
        // 1-on classe Pop_rec dans l'ordre des nouvelles dates de naissance:
        Collections.sort(Pop, new MonComparator());

        // 2- on initialise les données environementale
        //day_init_temp = -1;



        for (int p = 0; p
                < Pop.size(); p++) {
            Poissons Po = (Poissons) Pop.get(p);
            /* INUTILE si on enregistre pas l'env initiale.
            if (Po.day_init > day_init_temp && (Po.day_init % Dataset_EVOL.dt_HyMo == 0)) {
            try {
            Dataset_EVOL.load_data(Po.day_init);
            } catch (IOException e) {
            System.out.println("probleme dans ponte_env_init : " + e.getMessage());
            }
            }
             */
            //System.out.println("init()...");
            Po.init(Po.day_init);
            //  day_init_temp = Po.day_init;


        }
    }

// ------------------------------------------------------------------------------------------------------
    static boolean dans_la_zone(int i, int j) {
        return ((i > i_min) && (i < i_max) && (j > j_min) && (j < j_max));


    }

    static boolean pas_dans_la_zone_tampon(int i, int j) {
//        return ((i > (i_min - sponge_grid)) && (i < (i_max + sponge_grid)) && (j > (j_min - sponge_grid)) && (j < (j_max + sponge_grid)));
        return ((i > (i_min + sponge_grid)) && (i < (i_min + sponge_grid)) && (j > (j_min + sponge_grid)) && (j < (j_max - sponge_grid)));


    } //--------------------------------------------------------------------------
    // Relation de densite dependance pour le calcul du nombre d'oeuf par individu :

    int densite_dependence_locale(int d, int i, int j) {
// NE MARCHE PAS

        int nb_oeuf;
        // nb_nouveaux_super_indiv : nombre d'oeuf pondu par l'individu Po
        // D_strates[d][i][j] : nombre d'oeufs déja pondu à cette location par d'autres individus

        // le nombre d'oeuf pondu par individu est égal à Pop_max/Pop_rec pour stabiliser (Simulation.nbponte_par_indiv)
        // MULTIPLIé par une fonction fn de la densité, de forme arctangeante, afin  de réduire l'avantage des zstrates à haute densité :
        // Pour voir la geule de cette fonction dans matlab :
        // >> x=0:0.1:10
        // >> plot(x, (atan(D_max/2 - x)/atan(D_max/2))/2 + 0.5); grid
        // ou x est la densité dans la strate (D_strates[d][i][j])

        // Fonction de densité dépendence : ARCTANGEANTE
        // en forme de arc tangeante, avec le point d'inflexion en D_max/2
        // et normalisée de sorte à valoir 1 en 0 et 0 en D_max.



        double p = (double) (1.f / 2.f);


        double inflexion = ((double) D_max_spatiotemp) * p;


        double fn = 0.5f + (Math.atan(inflexion - D_strates_spatiotemp[d][i][j]) / Math.atan(inflexion)) / 2.0f;

        /*
        for (int d = 0; d<D_max; d++){
        double fn = 0.5f + (Math.atan(inflexion - d) / Math.atan(inflexion))/2;
        System.out.println(fn);
         */
        // Fonction de densité dépendence : LINEAIRE
        // >> plot(x, (-x + 79)/79); grid
        //double fn = (D_max- D_strates[d][i][j])/D_max;

        nb_oeuf =
                (int) Math.floor(Simulation.nbponte_par_indiv * fn) + 1;


        return nb_oeuf;


    }

    public float arondi(double nombre_long) {
        float nombre_cour = (float) Math.round(nombre_long * 100) / 100;


        return (nombre_cour);


    }

    public void trouver_strat_ponte_potentielles(Poissons Po) {


        int dd, d, ii, jj, i, j, explo_temp;


        double prof_fond, rayon_exploHE_maille;



        int k = 0; // on ne regarde que les temp et salt de surface.



        double Chla_ij = 0;


        boolean env_ok;

        // ON CREE le tableau pour que ce poisson sauvegarde ses zones de ponte :
        Po.init_strates();

//            // Strate temporelle la plus proche de Po.day_init, selon la resolution temporelle "resolution_temporelle_scanEnv" :
//            d = Math.round(Po.day_init / resolution_temporelle_scanEnv);

        // Convertir le rayon d'exploration de ce poisson en nombre de mailles :

        // HOMING TEMPOREL ?------------------------------------------------


        if (Simulation.HomingTemporel) {
            explo_temp = (int) Math.round(Po.rayon_exploration_temporel / Simulation.resolution_temporelle_scanEnv);


        } else {
            explo_temp = (int) 180.0f / Simulation.resolution_temporelle_scanEnv;


        }

        liste_d = Simulation.creer_liste_d(explo_temp, Po.day_init);
        // -----------------------------------------------------------------

        // HOMING GEOGRAPHIOQUE ?-------------------------------------------


        if (Simulation.HomingGeographique) {
            rayon_exploHE_maille = Po.rayon_exploration_spatial / Population.mean_grid_size;


        } else {
            rayon_exploHE_maille = 9999999.0f;


        } // HOMING TEMPERATURE ?  -------------------------------------------
        if (!Simulation.HomingTemperature) {
            Po.tolerance_temperature = 9999999.0f;


        } // HOMING SALINITE ?  ----------------------------------------------
        if (!Simulation.HomingSalinite) {
            Po.tolerance_salinite = 9999999.0f;


        } // HOMING HBL ?  ----------------------------------------------
        if (!Simulation.HomingHBL) {
            Po.tolerance_HBL = 9999999.0f;


        } // HOMING BATYMETRIE ?  ----------------------------------------------
        if (!Simulation.HomingBathymetrie) {
            Po.tolerance_batymetrie = 9999999;


        }

// HOMING GEOGRAPHIOQUE :
        for (ii = Math.max(Math.round((float) (Po.x_init - rayon_exploHE_maille)), Math.round((float) (x_min)));
                ii
                < Math.min(Math.round((float) (Po.x_init + rayon_exploHE_maille)), Math.round((float) (x_max))); ii++) {
            // PASSAGE AU COORDONEE DANS LE DOMAINE CHOISIT (de dimension X) :
            i = ii - Math.round((float) x_min);
            //System.out.println("i = " + i);


            for (jj = Math.max(Math.round((float) (Po.y_init - rayon_exploHE_maille)), Math.round((float) (y_min)));
                    jj
                    < Math.min((int) Math.round((float) (Po.y_init + rayon_exploHE_maille)), Math.round((float) (y_max))); jj++) {
                // PASSAGE AU COORDONEE DANS L DOMAINE CHOISIT (de dimension Y) : ---
                j = jj - Math.round((float) y_min);
                //System.out.println(" j = " + j);

                // VERIFICATION 1 : ce point est-il dans l'eau?


                if (Dataset_EVOL.isInWater(ii, jj)) {//&& !Dataset_EVOL.isCloseToCost(ii, jj)) {
                    // HOMING SUR LA BATHYMETRIE :
                    prof_fond = Dataset_EVOL.getDepth(ii, jj, 0);


                    if (Po.bathy_init > (prof_fond - Math.abs(Po.tolerance_batymetrie * prof_fond / 100))
                            && Po.bathy_init < (prof_fond + Math.abs(Po.tolerance_batymetrie * prof_fond / 100))) {

                        // HOMING TEMPOREL :
                        // intervalle temporel :
                        for (dd = 0; dd
                                < liste_d.length; dd++) {

                            int jour = liste_d[dd];

                            d = (int) Math.floor(jour / Simulation.resolution_temporelle_scanEnv);

                            // VERIFICATION 2 : La concentration en Chla en ce point/date est elle superieure au seuil choisi (nourriture adultes)
                            //                  OU bien on ne definis pas la zone de pont en fonction de la Chla


                            int mois = (jour / 30) % 12; // ON CONSIDERE QUE LES MOIS ONT 30 JOURS...



                            if (Simulation.Ponte_Chla_SeaWiFS) {
                                Chla_ij = Simulation.Chla[mois][i][j];


                            }

                            if ((Chla_ij > Simulation.Seuil_ponte_Chla_SeaWiFS) || !Simulation.Ponte_Chla_SeaWiFS) {

                                // VERIFICATION 3 : La densite de ponte a cette date n'est elle pas deja trop importante?
                                if ((!Simulation.flag_Densite_max) || (D_strates_temporelle[d] < Population.D_max_temporel)) {
                                    //System.out.println(" d = " + d);

                                    // HOMING SUR Température & Salinité (& HBL):
                                    env_ok = false;
                                    env_ok =
                                            ((Po.temperature_natal > Simulation.T[d][i][j][k] - Po.tolerance_temperature && Po.temperature_natal < Simulation.T[d][i][j][k] + Po.tolerance_temperature) && // HOMING TEMPERATURE
                                            (Po.salinite_natal > Simulation.S[d][i][j][k] - Po.tolerance_salinite && Po.salinite_natal < Simulation.S[d][i][j][k] + Po.tolerance_salinite));// && // HOMING SALINITE
//                                                    (Po.profThermo_natal > Th[d][i][j] - Po.tolerance_HBL && Po.profThermo_natal < Th[d][i][j] + Po.tolerance_HBL));//  HOMING PROF THERMO



                                    if (Simulation.flag_Densite_max && env_ok) {
                                        env_ok = Population.D_strates_spatiotemp[d][i][j] < Population.D_max_spatiotemp;


                                    }

                                    if (env_ok) {
                                        //System.out.println("env_OK");
                                        // System.out.println("ii = " + ii + " ; jj = " + jj + " ; d = " + d + " ; k = " + k);
                                        Po.Nb_strates++;
                                        // ON AJOUTE CETTE STRATE DANS LA LISTE DES Env_OK  DE CE POISSON :
                                        Po.strates[d][i][j] = true;

                                        // Alternative EVOL-DEB (car pb de heap memoire sinon) :
                                        // ne stocker que les d,i, j des env_ok :
                                        //Po.zones_de ponte_potentielles[Po.Nb_strates] = {d,i,j};



                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    static public void effacer_strates() {
        for (int s = 0; s
                < Pop.size(); s++) {
            Poissons Po = (Poissons) Pop.get(s);
            Po.erase_strates();


        }
    }

    public void stepSerial_TEST_DEB() {
        System.out.println("TEST DEB : individus reçoivent T et f constant - pas d'advection");

        // CORRECTION BUG ADVECTION (19 Juin 2009) pour interpolation :
        nb_jours_ecoule_entre_roms_records = -1;



        int jour_init = 0;


        int jour_final = Simulation.nb_jours_par_ans;

// -----------------------  PARCOUR DES JOURS DE L'ANNEE  ----------------------


        for (int jour = jour_init; jour
                < jour_final; jour++) {

            //System.out.println("nb_jours_ecoule_entre_roms_records = " +nb_jours_ecoule_entre_roms_records);
            nb_jours_ecoule_entre_roms_records++;


            if ((jour % 10) == 0) {
                System.out.println("TEST_DEB - GEN " + Simulation.compteur_year + " :  jour de l'annee: " + jour + " -------------- Pop.size = " + Pop.size());


            } // nb_jours_ecoule_entre_roms_records : Aconserver pour avoir le même pas de temps que dans le vrai run
            if (jour % Dataset_EVOL.dt_HyMo == 0) {
                nb_jours_ecoule_entre_roms_records = 0;


            }

            //**** ICI ON PARCOUR LA POPULATION --------------------------------
            for (int s = 0; s
                    < Pop.size(); s++) {
//                    System.out.println("Poisson " + ss + " début...");
                Poissons Po = (Poissons) Pop.get(s);



                if ((jour >= Po.day_init || Po.age > 0)
                        && Po.living) {
                    Po.stepSerial_TEST_DEB(jour); // <-- ADVECTION OU DEPLACEMENT


                } // PONTE ||
                if (Po.DEB.adult_ready_to_spawn) {
                    Po.DEB.spawn();


                }

            } // fin de la boucle sur les individus

            // ENREGISTREMENT DES POSITIONS
/*
            int sample_dt = 30; // jours
            if (jour % sample_dt == 0) {
            int mois = jour / sample_dt;
            // System.out.println("ENREGISTREMENT DES POSITIONS DU MOIS " + mois + " ... ");
            openPOSITION(Simulation.Strategie, mois);
            writePOSITION();
            closePOSITION();
            }
             *
             */

            if (Simulation.outpout_netcdf_flag) {
                //System.out.println(" REMPLISSAGE DU NetCDF TEST  " );
                try {
                    Sim_writer.writeResults_NetCDF(Population.Pop);


                } catch (Exception ex) {
                    System.out.println("probleme dans creation netcdf TES_DEB pop : " + ex.getMessage());


                }
            }
            Simulation.compteur_t_rec_intra_year++;



        } // fin de la boucle sur les jours
    }

    public void ponte_init_TEST_DEB() {
        // PONTE SEULEMENT LE JOUR 1
        int jour_debut_ponte = 0;


        int jour_fin_ponte = 1;



        for (int jour = jour_debut_ponte; jour
                < jour_fin_ponte;
                jour++) {
            //System.out.println("Simulation.nb_jours_par_ans = " + Simulation.nb_jours_par_ans);

            for (int p = 0; p
                    < Simulation.nbFishesJour; p++) { // Ici "Simulation.nbFishes_pasdetemps" a été remplacé par "Simulation.nbFishesJour" le 21 mars 2013
                Poissons Po = new Poissons(jour);
//                compteur_id_poissons++;
                Pop.add(Po);
            }
        }
    }



void RecMax(){

// superficie du plateau continental de 0 a 200m (config PEDRO)
float surf_plateau_nbmaille = 2199; // [en nombre de points de grille](environ)
// ATTENTION :  A RECALCULER SI ON CHANGE DE GRILLE
nrec_max_jour = 1;
// on veut un recrument max de nrec_max_jour individu par jour
// donc il faut un recmax de (nrec_max_jour/surf_plateau_nbmaille)*surface_zone dans chaque zone


// Nombre de zones (de 1 degre de latitude) :
 nbzones = (int) Math.round(Dataset_EVOL.lat_max - Dataset_EVOL.lat_min);

// soit N_recmax_SI le nombre max de recrue par zone  :
N_recmax_SI = new float[nbzones];

int[][] mask_plateau;
int[][] mask_zone;
int[][] mask_latmin_zone, mask_latmax_zone;
mask_plateau = new int[Dataset_EVOL.ny][Dataset_EVOL.nx];
mask_zone = new int[Dataset_EVOL.ny][Dataset_EVOL.nx];
mask_latmin_zone = new int[Dataset_EVOL.ny][Dataset_EVOL.nx];
mask_latmax_zone = new int[Dataset_EVOL.ny][Dataset_EVOL.nx];

        for (int j = 0; j< Dataset_EVOL.ny; j++) {
            for (int i = 0; i < Dataset_EVOL.nx; i++) {
                mask_plateau[j][i] = (Dataset_EVOL.hRho[j][i] < Dataset_EVOL.prof_talu) ? 1:0;
                mask_plateau[j][i] = mask_plateau[j][i] * Dataset_EVOL.maskRho[j][i];
            }
        }


for(int zone=0; zone<nbzones; zone++){
    // Surface en nombre de maille entre lat_min et lat_min+1 :
    float surf_grd_zone_i = 0;

        for (int j = 0; j< Dataset_EVOL.ny; j++) {
            for (int i = 0; i < Dataset_EVOL.nx; i++) {
                mask_latmin_zone[j][i] = (Dataset_EVOL.latRho[j][i] > (Dataset_EVOL.lat_min + zone)) ? 1:0;
                mask_latmax_zone[j][i] = (Dataset_EVOL.latRho[j][i] < (Dataset_EVOL.lat_min + zone+1)) ? 1:0;
                mask_zone[j][i] = mask_plateau[j][i] * mask_latmin_zone[j][i] * mask_latmax_zone[j][i];

                surf_grd_zone_i += mask_zone[j][i];
            }
        }
    
//    N_recmax_SI[zone]  = (int) (Math.floor(surf_grd_zone_i*(nrec_max_jour/((float) surf_plateau_nbmaille))));
    N_recmax_SI[zone]  = surf_grd_zone_i*(nrec_max_jour/((float) surf_plateau_nbmaille));

   System.out.println("Zone " + zone + "  --> N_recmax_SI = " + N_recmax_SI[zone]);
}
} // end of class
}

