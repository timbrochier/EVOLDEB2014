package evolution;

import java.io.*;
import java.util.*;

import util.GetChla_SeaWiFSclim;
import util.GetOxyclinDepth_BIDOUILLE;
import util.Get_SpawnFile;
import util.GetChlaPISCES_integ;
import util.GetCarbonZPISCES_integ;
import util.Lire_configuration;
import util.OutpoutManager_netcdf;

public class Simulation extends Object {
// --------- DECLARATION DES VARIABLES -----------------------------------------

    public static int compteur_year, dt_advec, compteur_t_rec_intra_year,Pop_max,nb_jours_par_ans;
    static int nbFishesJour, year, ageMinAtRecruitment, nb_strat_temporelles, Nb_strates_total, nbFishes_pasdetemps;
    public static float temp_lethal_min, temp_lethal_max, patchiness_rad_HS, vertical_patchiness_rad, patchiness_rad_OPP, nbponte_par_indiv;
    // Range de tolerances :
    public static float rayon_exploration_spatial_max, rayon_exploration_spatial_min, rayon_exploration_temporel_max, rayon_exploration_temporel_min, Range_tolerance_temperature, Range_tolerance_salinite, Range_tolerance_batymetrie, Range_tolerance_HBL;
    static int time_rad_HS, nbiter, resolution_temporelle_scanEnv, Nb_oeufs_pondu;// alea_temps
    // Range des parametres environnnementaux (lors du scan de lenvironnement):
    static float[] Temp_min, Temp_max, Salt_min, Salt_max;
    public static int egg_duration, nb_record;
    static double facteur_ponte_strate;
    //double dist_natal, rayon_exploHE;
    FileWriter resultsWriter, paramWriter;
    Dataset_EVOL dataset_EVOL;
    OutpoutManager_netcdf netcdf_outpout_writer;

    Population population;
    static double seuil_dispersion, dt, egg_excess_buoyancy, superindiv_mort_moyen, decalage_zeta, seuil_concentration, Swim_lowerlimit, Swim_upperlimit;
    public static double Seuil_ponte_Chla_SeaWiFS, Seuil_rec_Chla_SeaWiFS, Seuil_rec_Chla_PISCES, Seuil_rec_CarbZ_PISCES;
    public static boolean flagRET_Chla_PISCES, flagGREG_dens, flagRET_shelf, flagRET_Chla_SeaWiFS, BetHedging, flagSwim_O2, flagSwim_random, record_initial_env, flagBuoy, flagGrowth, flagTurbulence, lire_HBL, flag_Densite_max, Ponte_Chla_SeaWiFS;
    public static boolean flagRET_CarbZ_PISCES, outpout_netcdf_flag, TEST_KINESIS;
    //public static int[] Nb_strates_poisson; //Nb_strates_date
    static int nbre_generations;
    public static String Strategie, output_dir;
    public static String dist_unit;
    static double nbeggs_max_per_SI ;
    boolean year_test_UNIQUE;
    static boolean ADVECTION;
        // pour tester DEB tout seul, avec T et f donnés dans poisson...
    boolean TEST_DEB;

    // OPPORTUNISME
    // cible
    static double Topp, Sopp, bathymax_opp, bathymin_opp, Th_opp, DensOpp, F_annuel, Min_fishing_size;
    // marge (=tolerances)
    static double dens_marge, salt_marge, temp_marge, Th_marge;
    static boolean[][][][] strates_opp;
    public static float[][][][] T; // temperature
    public static float[][][][] S; // salinidad
    public static float[][][] Th;  // prof de thermocline
    public static double[][][] Chla;  // Chlorophylle de surface (SeaWiFS)
    int[] liste_year_TEST_CRIT_REC;
    public static boolean HomingGeographique, HomingTemporel, HomingTemperature, HomingSalinite, HomingBathymetrie, HomingHBL;
    
    // --------------------- END OF VARIABLES DECARATIONS ----------------------
    public Simulation() {
        //-------------- PARAMETRES A DEFINIR PAR L'UTILISATEUR ------------------------

// A -  LA REGION
        dataset_EVOL = new Senegal(); //<-- ICI METTRE Peru OU CANARY OU Chile

outpout_netcdf_flag = true;
if(outpout_netcdf_flag){
        netcdf_outpout_writer = new OutpoutManager_netcdf();
}
// DANS EVOL 2010,  la pluspart des parametres sont lu dans le fichier config :

// PARAMETRES LU DANS LE FICHIER CONFIG
// [EVOL2010] Lors du lancement d'une simulation, donner en argument un fichier
// config (ascii) contenant la liste des parametres et leur valeur
// ex : Strategie = BH
//      nbre_generations=30
//      etc..
// Donner ce fichier config en argument lors de l'execution :
// java -Xmx4g EVOL2010.jar fichier_config.txt


// Pour une ponte selon les observations (complete par climato), activer la clef
// "PONTE_OBS" dans la presente classe (Simulation).

        try {

// B - STRATEGIE DE REPRODUCTION :
//        Strategie = "BH"; // HS, HE, Opp, BH
            // HS = Homing Strict (ou geographique)
            // HE = Homing Environnemental
            // Opp= Opportunisme
            // BH = Bet Hedging
            Strategie = Dataset_EVOL.param_config.getConfig_string("Strategie");
//    SIMU= (int) param_config.getConfig_int("SIMU");
            System.out.println("Strategie = " + Strategie);


// C - NOMBRE DE GENERATIONS
            nbre_generations = Dataset_EVOL.param_config.getConfig_int("nbre_generations"); // Nombre d'annees ou se répete la strategie de reproduction

// D - NOMBRE DE POISSONS INITIALEMENT PONDU PAR JOUR
            nbFishesJour = Dataset_EVOL.param_config.getConfig_int("nbFishesJour");

// E - CRITERES PHYSIQUES AFECTANT LE RECRUTEMENT
            ageMinAtRecruitment = Dataset_EVOL.param_config.getConfig_int("ageMinAtRecruitment");   // [jours]  Age des individus pour les test de recrutement (= duree de l'advection des larves)

            // ZONE DE RETENTION -------------
            flagRET_shelf = Dataset_EVOL.param_config.getConfig_boolean("flagRET_shelf");            // Recrutement si position finale est sur le shelf

            // DENSITE DEPENDANCE -------------
            flag_Densite_max = Dataset_EVOL.param_config.getConfig_boolean("flag_Densite_max");          // Pour fixer une densité maximum d'individu par strate lors de la ponte (voir Population.zoneGeo2grid)

// D - PARAMETRES DE LA STRATEGIE DE PONTE

            // Patchiness
            vertical_patchiness_rad = Dataset_EVOL.param_config.getConfig_int("vertical_patchiness_rad");// [metres]         ALL STRATEGY (prof de ponte = prof de naissance +- verticalpatch_HS)
            patchiness_rad_OPP = Dataset_EVOL.param_config.getConfig_int("patchiness_rad_OPP");     // [Kilometres]        OPPPORTUNISM

            resolution_temporelle_scanEnv = Dataset_EVOL.param_config.getConfig_int("resolution_temporelle_scanEnv");            // [jours] pour la fonction scan_env ; on scan l'environement tout les x jours
            //alea_temps = 10;        // [jours] les adultes rechercheent leur environement de naisance à leur date de naissance + ou - dt_jour jours

            // TOLERANCES
            // Tolerence autour de l'environnement recherche autour de celui de naissance (+ ou - cette valeur)
            // chaque individu se differencie par ses tolerences (+ ou -), tiree aleatoire au début de la simu entre les valeurs min et max choisies
            // Spatial : Rayon de patchiness horizontal en Km (depuis le 24 Juin 2009, avant : en mailles)

            rayon_exploration_spatial_min = Dataset_EVOL.param_config.getConfig_int("rayon_exploration_spatial_min"); // [Kilometres]
            rayon_exploration_spatial_max = Dataset_EVOL.param_config.getConfig_int("rayon_exploration_spatial_max"); // [Kilometres]

            // Temporel : les individus pondent à leur date de naissance +- time_rad_HS jours
            rayon_exploration_temporel_min = Dataset_EVOL.param_config.getConfig_int("rayon_exploration_temporel_min");  // [Nb Jours]
            rayon_exploration_temporel_max = Dataset_EVOL.param_config.getConfig_int("rayon_exploration_temporel_max");  // [Nb Jours]
            Range_tolerance_temperature = (float) Dataset_EVOL.param_config.getConfig_double("Range_tolerance_temperature");//0.2f;// [°C] Max de tolerance en temperature (min =0):
            Range_tolerance_salinite = (float) Dataset_EVOL.param_config.getConfig_double("Range_tolerance_salinite");//0.04f; // [PSU] Max de tolerance en Salinite (min =0):
            Range_tolerance_batymetrie = Dataset_EVOL.param_config.getConfig_int("Range_tolerance_batymetrie");        // [% de la profondeur "home"]
            Range_tolerance_HBL = Dataset_EVOL.param_config.getConfig_int("Range_tolerance_HBL");               // [metres] NON UTILISEE (aout 2009)

            // HOMING ENVIRONNEMENTAL : QUELLE VARIABLES POUR LE HOMING ?
            HomingGeographique = Dataset_EVOL.param_config.getConfig_boolean("HomingGeographique");
            HomingTemporel = Dataset_EVOL.param_config.getConfig_boolean("HomingTemporel");
            HomingTemperature = Dataset_EVOL.param_config.getConfig_boolean("HomingTemperature");
            HomingSalinite = Dataset_EVOL.param_config.getConfig_boolean("HomingSalinite");
            HomingBathymetrie = Dataset_EVOL.param_config.getConfig_boolean("HomingBathymetrie");
            HomingHBL = Dataset_EVOL.param_config.getConfig_boolean("HomingHBL"); // Pas pertinent car Couche de mélange varie pour des raisons differentes selon les zones (upwelling, stratification,...)


            // OPPORTUNISME
            // 1 - environnement cible
            DensOpp = Dataset_EVOL.param_config.getConfig_double("DensOpp"); // Densite de l'eau visee par opportunistes
            Sopp = Dataset_EVOL.param_config.getConfig_double("Sopp"); // Salinite de l'eau visee par opportunistes
            Topp = Dataset_EVOL.param_config.getConfig_double("Topp"); // Temperature de l'eau visee par opportunistes
            Th_opp = Dataset_EVOL.param_config.getConfig_int("Th_opp"); // prof thermocline  visee par opportunistes
            bathymax_opp = Dataset_EVOL.param_config.getConfig_int("bathymax_opp"); // Bathymax opportunistes
            bathymin_opp = Dataset_EVOL.param_config.getConfig_int("bathymin_opp"); // Bathymin  visee par opportunistes

            // 2 - tolerances de variation autour de l'environnement cible
            dens_marge = Dataset_EVOL.param_config.getConfig_double("dens_marge"); // TOLERANCE de Densite de l'eau visee par opportunistes
            salt_marge = Dataset_EVOL.param_config.getConfig_double("salt_marge"); // TOLERANCE de Salinite de l'eau visee par opportunistes
            temp_marge = Dataset_EVOL.param_config.getConfig_double("temp_marge"); // TOLERANCE de Temperature de l'eau visee par opportunistes
            Th_marge = Dataset_EVOL.param_config.getConfig_int("Th_marge"); //TOLERANCE de prof thermocline  visee par opportunistes


            // E - PROPRIETE DES INDIVIDUS
            temp_lethal_min = (float) Dataset_EVOL.param_config.getConfig_double("temp_lethal_min");      // [°C]  temperature lethale minimum
            temp_lethal_max = (float) Dataset_EVOL.param_config.getConfig_double("temp_lethal_max");      // [°C]  temperature lethale maximum
            flagGrowth = Dataset_EVOL.param_config.getConfig_boolean("flagGrowth");           // (Growth ultra simple de ichthyop (Koné, 2006) A FAIRE : adapter modele de Urtizberea et al. 2009
            flagBuoy = Dataset_EVOL.param_config.getConfig_boolean("flagBuoy");              // [boolean] flotabilité des oeuf: si false neutre, si true
            egg_excess_buoyancy = Dataset_EVOL.param_config.getConfig_double("egg_excess_buoyancy"); // [kg.L-1] // d'après Gorant et al. 2007
            egg_duration = Dataset_EVOL.param_config.getConfig_int("egg_duration"); // [jours] durée du stade oeuf
            flagSwim_O2 = Dataset_EVOL.param_config.getConfig_boolean("flagSwim_O2");             // DVM entre surface et oxycline... a valider avec J. Tam
            flagSwim_random = Dataset_EVOL.param_config.getConfig_boolean("flagSwim_random");
            Swim_upperlimit = Dataset_EVOL.param_config.getConfig_double("Swim_upperlimit");
            Swim_lowerlimit = Dataset_EVOL.param_config.getConfig_double("Swim_lowerlimit");

            
            // F - PECHE
         F_annuel = (float) Dataset_EVOL.param_config.getConfig_double("F_annuel");      // [taux] % de poissons pechés par ans
         Min_fishing_size = (float) Dataset_EVOL.param_config.getConfig_double("Min_fishing_size");      // [cm] taill min pour les poissons pêchés
         
         
    
            
// G - DIVERS
            dist_unit = "Km";             // Unite pour les distances de patchiness (Km = kiometres ou Nm = Miles nautiques)
            flagTurbulence = Dataset_EVOL.param_config.getConfig_boolean("flagTurbulence");       // Pour prendre en compte la turbulence sous-maille dans le transport lagrangien
            lire_HBL = Dataset_EVOL.param_config.getConfig_boolean("lire_HBL");             // Pour lire la couche de mélange HBL (nécéssaire pour Homing Env si on veut)
            nb_jours_par_ans = Dataset_EVOL.param_config.getConfig_int("nb_jours_par_ans");       // Nombre de jours par ans dans la simu ROMS (360 dans certaines simus, 365 dans les nouvelles)
            dt_advec = Lire_configuration.getConfig_int("dt_advec");                 // [minutes] pas de temps pour l'advection des larves (interpolation des champs de vitesse tout les dt_advec heures)
            // IMPORTANT : choisir dt_advec tel qu'on ait un nombre entier de dt_advec entre deux sorties roms)
            record_initial_env = Dataset_EVOL.param_config.getConfig_boolean("record_initial_env");   // Pour enregistrer les T, S et HBL du lieu de naissance (augmente temps de simulation). Automatique pour le Homing environemental et l'Opportunisme
            nbeggs_max_per_SI = 1e14;


        } catch (Exception e) {
// TODO Auto-generated catch block
            System.out.println("probleme dans Lire_configuration : " + e.getMessage());
            e.printStackTrace();
        }

        /*
        if (Strategie.equals("HE") || Strategie.equals("Opp") || flagBuoy) {
        record_initial_env = true;
        }
         */
        //... et parceque ça simplifie VRAIMENT (mais vraiment) TOUT (par rapport a 365 jours) :
        if (resolution_temporelle_scanEnv > 1) {
            nb_jours_par_ans = 360;
        }

        TEST_DEB= false; // true : T et bouffe constant, pas de spatial
        TEST_KINESIS = false; // true : champs physique et bioch constants
        ADVECTION = true; // false : enlève l'advection des courants 3D pour les adultes

        //           if (Strategie.equals("HS")) {
        //       resolution_temporelle_scanEnv = 1;
        //   }
//--------------- END OF USER'S CHOICES -----------------------------------
        decalage_zeta = 0; // Pour recaler les profondeur quand on a un décalage (ex : run preindustrielle)
// PARTICULARITE DE CERTAINES CONFIG :

        if (Dataset_EVOL.region.equals("PERU_PISCES_sixieme") || Dataset_EVOL.region.equals("PERU_sixieme")) {
            nb_jours_par_ans = 365;
        } else if (Dataset_EVOL.region.equals("PERU_preindus_rv10")) {
            nb_jours_par_ans = 360;
            decalage_zeta = 6.566;
            System.out.println("decalage ZETA pour = 6,74m");
        } else if (Dataset_EVOL.region.equals("PERU_4xCO2")) {
            nb_jours_par_ans = 360;
        }
    }

    public void initSerial() {

        // TEST creer_liste_d
/*        liste_d = creer_liste_d(3, 330);
        for (i = 0; i<liste_d.length; i++){
        System.exit(0);
        // FIN TEST
         */

        System.out.println("EVOL 2014 : init serial ... simu " + Dataset_EVOL.SIMU);
        System.out.println("Config file =  " + IBM.config_file);

        System.out.println("Peche : F_annuel = " + F_annuel + " ; soit F = " + F_annuel / (float) nb_jours_par_ans + " / jour ");
        
        
        nbFishes_pasdetemps = nbFishesJour * resolution_temporelle_scanEnv;

        nb_strat_temporelles = nb_jours_par_ans / resolution_temporelle_scanEnv;

        Pop_max = nbFishesJour * nb_jours_par_ans;


if (TEST_KINESIS){
    Pop_max = nbFishesJour;
}
        // Initialisation du compteur de generation a 0 (1 generation = 1 year)
        compteur_year = 0;  // <--------------------- POUR ALLER PLUS VITE EN FIN DE THESE, LE REMETRE A ZERO APRÈS!!!!!!!!!!!!!

        // Selection de la 1er année hydro dans la liste aléatoire
        year = Dataset_EVOL.yearlist100[compteur_year];


        // Lecture des champs constants dans les sorties roms (dans le fichier mois 1 de l'année "year" selectionnee ci-dessus)
        try {
            dataset_EVOL.setUp();
        } catch (IOException e) {
            System.out.println("Probleme dans setup : " + e.getMessage());
        }

        if (Ponte_Chla_SeaWiFS || flagRET_Chla_SeaWiFS) {
            try {
                GetChla_SeaWiFSclim.setUp();
            } catch (IOException e) {
                System.out.println("Probleme dans GetChla_SeaWiFSclim.SetUp : " + e.getMessage());
            }
        }

        if (flagRET_Chla_PISCES) {
            try {
                GetChlaPISCES_integ.setUp();
            } catch (IOException e) {
                System.out.println("Probleme dans GetChla_SeaWiFSclim.SetUp : " + e.getMessage());
            }
        }

        if (flagRET_CarbZ_PISCES) {
            try {
                GetCarbonZPISCES_integ.setUp();
            } catch (IOException e) {
                System.out.println("Probleme dans GetCarbonZPISCES_integ.SetUp : " + e.getMessage());
            }
        }
        if (flagSwim_O2) {
            try {
                GetOxyclinDepth_BIDOUILLE.setUp();
            } catch (IOException e) {
                System.out.println("Probleme dans GetOxyclinDepth_clim.SetUp : " + e.getMessage());
            }
        }

        if (flagBuoy) {
            util.BuoyancyScheme.init();
        }

        if (flagGrowth) {
            util.GrowthModel.init();
        }

        // Create the output directory
        output_dir = (Dataset_EVOL.main_output_directory + Dataset_EVOL.region + "/" + Strategie + "_" + Dataset_EVOL.region + "_simu" + Dataset_EVOL.SIMU + "/");
        boolean success = (new File(output_dir).mkdir());
        if (success) {
            System.out.println("Directory: " + output_dir + " created");
        }
        System.out.println("Enregistrement des sorties dans " + output_dir);

        // Creation de la population
        System.out.println("new population ...");
        population = new Population(netcdf_outpout_writer);

        // Calcul de la zone de ponte initiale dans la grille (transfer de lat, lon en i, j)
        population.ZoneGeo2Grid();
        
        // Calcul du nombre de recrue max par degre de latitude de cote
        population.RecMax();

        // TEST DEB debut

        if (TEST_DEB){
TEST_DEB();
          }

// TEST DEB fin

        System.out.println("Ecriture des parametres de strategie dans le fichier param.txt");
        openParam();
        closeParam();

        // Lancement de la ponte initiale
        System.out.println("Nombre de generations : " + nbre_generations);


PonteInitiale();

    }

//---------------------------------------------------------------------------------------------------------------
    public void PonteInitiale() {

        // 1er lecture des champs dynamique
        try {
            dataset_EVOL.init();
        } catch (IOException e) {
            System.out.println("blem 1: " + e.getMessage());
            System.exit(0);
        }


        // Calcul du nombre d'iteration a faire lors de l'advections des individus (OEUFS ET LARVES):
        nbiter = (int) (((double) Dataset_EVOL.dt_HyMo * 24 *60) / (double) dt_advec); //(dt_HyMo est en jour et dt_advec en MINUTES)
        
        System.out.println("Pas de temps pour l'advection des larves : " + dt_advec/60 + "heures");
        System.out.println("Soit " + nbiter + " interpolations entre les sorties roms.");
        // nbre d'iteration par pas de temps des sorties roms         // exemple : peru : sorties a 2 jours--> 12 iteration = une actualisation de champs de vitesse toute les 4 heure
        System.out.println("--------------------->  PONTE INITIALE  <--------------------");
        System.out.println("Nbre d'individu par jour : " + nbFishesJour);

            population.ponte_init();

        System.out.println("    xxxxxxxxxxxxxxxxxxx  -  ponte initiale terminee   - xxxxxxxxxxxxxxxx");

        if (Strategie.equals("HE")){
        System.out.println("Scan de l'environement: on lit les conditions de Température, Salinité, et profondeur de thermocline pour chaque point de grille...");
        //System.out.println("Les profondeurs scannées sont : " + Dataset_EVOL.prof_potentielles);
        scan_env();
        System.out.println("Scan de l'environement finis.");
        }

        // CREATION DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" CREATION DU NetCDF TEST  " );
            try {
            netcdf_outpout_writer.openResults_NetCDF();
            netcdf_outpout_writer.open_Biom2D_NetCDF();
            netcdf_outpout_writer.open_DataIndiv_NetCDF();          
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf 2 : " + ex.getMessage());
            }
}
        // ---------------------------------------------------------------------

        // A ce point, tout les poissons initiaux de la premiere annee sont en place
        // maintenant on les fait vivre et mourrir :

        System.out.println("--------------@@@@ DEPLACEMENT INITIAL  @@@@------------------");
        compteur_t_rec_intra_year =0;
        population.stepSerial();
        System.out.println("-----******------ FIN DU DEPLACEMENT INITIAL ------******-----");


        //**********************************************

        // Test de recrutement (pour les contraintes qui se vérifient à la fin , GREG et RET)
        //population.test_recruit();

        // Enregistrement des résultats :
        System.out.println(" ----> Ecriture des zones des positions de pontes, et finales ...");
        
        // FERMETURE DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" FERMETURE DE L'ANCIEN  NetCDF  " );
            try {
            netcdf_outpout_writer.CloseResults_NetCDF();
            netcdf_outpout_writer.Close_Biom2D_NetCDF();
            netcdf_outpout_writer.Close_DataIndiv_NetCDF();           
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf 2 : " + ex.getMessage());
            }
}
        // ---------------------------------------------------------------------

        //ancienne sorties (ASCII)
        //openResults(Strategie);
        //writeResults();
        //closeResults();

        System.out.println(" Ecriture finie.");

        // On copie la liste des individu recruté dans la liste "Pop_rec":
//        population.copieListe(); <-- (EFFACE LA LISTE "Pop")

        // On remet a zero les compteurs de densité :
        if (flag_Densite_max || flagGREG_dens) {
            population.mise_a_zero_densites();
        }


        // ANNEE SUIVANTE :
        // On reprend pareil sauf qu'il n'ya plus de distribution initiale
        // de la ponte, on aurra des juveniles et adultes qui font de la
        // nouvelle ponte en fonction de leur energie et strategie :

        // Homing Strict : ponte de patche autour des recrutes (x_init, y_init, z_init)
        // Homing Envir : chaque recrute pont au hazard avec outzone = temp> temp_init ou salt > salt_init, etc...
        // Opportunism : chaque recrute pont au hazard avec outzone = temp> temp_opp ou Salt> Salt_opp , etc...
        // Bet Hedging : chaque recrute pond partout au hazard dans la zone initiale? ou bien  comme HS mais avec un rayon de patche tres grand (BH spatial)

        // Lancer la recurssion :
        compteur_year = compteur_year + 1; // depuis le 4 octobre 2010
        Boucle();
    }
// -----------------------------------------------------------------------------

    public void Boucle() {

        while (compteur_year < nbre_generations) {
            year = Dataset_EVOL.yearlist100[compteur_year];

            System.out.println("<<<<<<<<<<<<<<<< GENERATION " + compteur_year + ">>>>>>>>>>>>>>>>");
            System.out.println("Annee de forcage physique : " + year);

// STRATEGIE DE PONTE :  //
            if (Strategie.equals("HE")) {
                HomingEnvironmental();
            }

            /*
             * NE FAIRE LA RECHERCHE DES STRATES ENVIR(HOMING) QUE POUR LES INDIV MATURE  ??
            if (PONTE_OBS) {
            SpawningObs();
            } else if (Strategie.equals("HS")) {
            HomingStrict();
            } else if (Strategie.equals("HE")) {
            HomingEnvironmental();
            } else if (Strategie.equals("Opp")) {
            Opportunism();
            } else if (Strategie.equals("BH")) {
            BetHedging();
            } else {
            System.out.println("Faut preciser la STRATEGIE!!!");
            System.exit(0);
            }

             */
        // CREATION DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" CREATION DU NetCDF TEST  " );
            try {
            netcdf_outpout_writer.openResults_NetCDF();
            netcdf_outpout_writer.open_Biom2D_NetCDF();
            netcdf_outpout_writer.open_DataIndiv_NetCDF();          
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf 3 : " + ex.getMessage());
            }
}
        // ---------------------------------------------------------------------

/////////////////////// DEPLACEMENT /////////////////////////////
            System.out.println("GEN " + compteur_year + "  -------******* Vie larvaire puis juvenile puis adulte... (deplacements)");
            compteur_t_rec_intra_year =0;
            population.stepSerial();
            System.out.println("GEN " + compteur_year + "  --------xxxx Fin de l'annee " + year + "..");

// Teste de recrutement (pour les contraintes qui se vérifient à la fin , GREG et RET)
//            population.test_recruit();

// ENREGISTREMENT DES ZONES DE PONTES ET POSITIONS FINALES            
//-->!! à modifier pour n'enregistrer que les ponte de cette année (faire une liste?)
            //        openResults(Strategie);
        // FERMETURE DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" FERMETURE DU NetCDF TEST  " );
            try {
            netcdf_outpout_writer.CloseResults_NetCDF();
            netcdf_outpout_writer.Close_Biom2D_NetCDF();
            netcdf_outpout_writer.Close_DataIndiv_NetCDF();            
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf 4 : " + ex.getMessage());
            }
}
        // ---------------------------------------------------------------------

// COPIE DE LA LISTE DES RECRUTES dans la liste Pop_rec
//            population.copieListe();
            System.out.println(" ------------ Generation " + compteur_year + "terminee --------------");
//
// On remet a zero les compteurs de densité :
            if (flag_Densite_max || flagGREG_dens) {
                population.mise_a_zero_densites();
            }
// On remet a zero les zones de ponte potentielles fde l'année:
            Population.effacer_strates();

// DEBUG : stat de sorties des adultes
            System.out.println(" sortie_NORD : " + Population.sortie_NORD);
            System.out.println(" sortie_SUD : " + Population.sortie_SUD);
            System.out.println(" sortie_EST : " + Population.sortie_EST);
            System.out.println(" sortie_OUEST : " + Population.sortie_OUEST);

            // ON REMET A ZERO : 
        Population.sortie_NORD = 0;
        Population.sortie_SUD = 0;
        Population.sortie_EST = 0;
        Population.sortie_OUEST = 0;

            compteur_year = compteur_year + 1;
        } // Fin du while (compteur_year < nbre_generations)
        if (compteur_year == nbre_generations) { // en principe ça devrait être le cas!!
            System.out.println(" c'est fini pour : " + Strategie);
            System.exit(0);
        }
    }

//========================================================================================================   @

//--------------------------------------  ++  HOMING STRICT ++ ----------------------------------------- :   @
    public void HomingStrict() {
        System.out.println("------%%%%%%----*****---- HOMING STRICT ----******-------%%%%%%%-------");

        // Calcul du nombre de ponte moyen par individu cette annee :
        nbponte_par_indiv = (float) Pop_max / Population.Pop_rec.size();

        // Et dans le cas des super individus : c'est variable pour chaque individu
        // Nbre moyen  de super_individumort :
//        superindiv_mort_moyen = Population.superindiv_mort_total / Population.Pop_rec.size();
//        System.out.println(" superindiv_mort_moyen : " + superindiv_mort_moyen);

        System.out.println(" nombre d'adultes pondeurs cette annee : " + Population.Pop_rec.size() + " , soit une moyenne de " + nbponte_par_indiv + " oeuf par individu");
        System.out.println("Ponte sur les traces des survivants de l'annee precedante...");

        // Ponte des individus a leur lieu et date de naissance, +- patchiness_rad_HS et time_rad_HS
        population.ponte_HomingStrict();
    }

//--------------------------------------- ++ HOMING ENVIRONEMENTAL  ++ --------------------------------- :   @
    public void HomingEnvironmental() {
        System.out.println("-----%%%%%-----*****----- HOMING ENVIRONMENTAL -----*****------%%%%%------");

        System.out.println("Dimensions du tablea de zones de ponte : ");
        System.out.println("nb_strat_temporelles = " + nb_strat_temporelles);
        System.out.println("X = " + Population.X + " , Y = " + Population.Y);
        System.out.println("prof_potentielles.length = " + Dataset_EVOL.prof_potentielles.length);

        //System.out.println("Population.Pop_rec.size() = " + Population.Pop_rec.size());
        // Calcul du nombre de ponte moyen par individu cette annee :
        //nbponte_par_indiv = (float) (Simulation.nbFishesJour * nb_jours_par_ans) / Population.Pop_rec.size();
        nbponte_par_indiv = (float) (Simulation.nbFishesJour * nb_jours_par_ans) / Population.Pop_rec.size();
        System.out.println("nbponte_par_indiv " + nbponte_par_indiv);

        System.out.println("Scan de l'environement: on lit les conditions de Température, Salinité, et profondeur de thermocline pour chaque point de grille...");
        //System.out.println("Les profondeurs scannées sont : " + Dataset_EVOL.prof_potentielles);

        scan_env();

        System.out.println("Scan de l'environement finis.");

        System.out.println("Recherche des strates environentales pour les quelles on a eu des survivants l'annee precedantes, parmis toute les strates potentielle....");

        int p, dd, d, ii, jj, i, j, k, explo_temp, Nb_indiv_ayant_pondu;
        float pourc;
        double prof_fond, rayon_exploHE_maille;
        double Chla_ij = 0;
        boolean env_ok;
        Nb_strates_total = 0;

        Nb_indiv_ayant_pondu = 0; // Compteur du nombre de poisson ayant pondu
        Nb_oeufs_pondu = 0;
    }


// Pour generer la liste des dates de ponte potentielles a partir de la date de naissance et du rayon d'exploration
   public static int[] creer_liste_d(int explo_temp, int jour_centre) {
        //System.out.println(" explo_temp  = " + explo_temp);
        //System.out.println(" jour_centre  = " + jour_centre);
        int date1, date2, i, j, k, d;

        date1 =
                jour_centre - explo_temp * resolution_temporelle_scanEnv;
        date2 =
                jour_centre + explo_temp * resolution_temporelle_scanEnv;
        int [] liste_d =
                new int[2 * explo_temp + 1];

        // Si la date1 se trouve AVANT le 1er janvier :
        if (date1 < 0) {
            date1 += nb_jours_par_ans;
            i =
                    0;
            for (d = date1; d
                    < nb_jours_par_ans; d =
                            d + resolution_temporelle_scanEnv) {
                liste_d[i] = d;
                i++;
            }


            for (d = 0; d
                    <= date2; d =
                            d + resolution_temporelle_scanEnv) {
                liste_d[i] = d;
                i++;

            }
            //OU BIEN Si la date2 se trouve APRES le 31 decembre :

        } else if (date2 >= nb_jours_par_ans) {
            date2 -= nb_jours_par_ans;
            i =
                    0;
            int last_d = 0;
            for (d = date1; d
                    < nb_jours_par_ans; d =
                            d + resolution_temporelle_scanEnv) {
                liste_d[i] = d;
                i++;

                last_d =
                        d;
            }
//1er jour en janvier :

            int first_jan = (last_d + resolution_temporelle_scanEnv) - nb_jours_par_ans;
            for (d = first_jan; d
                    <= date2; d =
                            d + resolution_temporelle_scanEnv) {
                liste_d[i] = d;
                i++;
            }
            // OU BIEN SI date1 et date2 au millieu de l'année :

        } else {
            i = 0;
            for (d = date1; d
                    <= date2; d =
                            d + resolution_temporelle_scanEnv) {
                liste_d[i] = d;
                i++;
            } 
        }
        return (liste_d);
    }

//--------------------------------------- ++ BET-HEDGING  ++ ------------------------------------------- :   @
    public void BetHedging() {
        System.out.println("-----%%%%%%-----****--- BET HEDGING ----********----------%%%%%%%%----------");

        // BetHedging = false; // pour pondre que dans la zone cotiere... // A VOIR
        System.out.println("Ponte dans toute la zone d'Upwelling");
        population.ponte_init();

    }
//========================================================================================================   @

//------------------------------------------------------------------------------------------------------ :   @
    public void scan_env() {
        System.out.println(" Temperature et salinites de chaque point de grille (strates)...");
        // convertir la zone de ponte potentielle en coordonnees grille : 

        // On pose: nombre de strates temporeles= nbre de sorties du modèles hydro (sur l'année):
        System.out.println("nb_strat_temp = " + nb_strat_temporelles);

        //	 Zones potentielle = zone initiale (x_min, y_min, la_min, lat_max))
        //	 initialisation des tableau de donnee environementale (T, S, Th):
        T =
                new float[nb_strat_temporelles][Population.X][Population.Y][Dataset_EVOL.prof_potentielles.length]; // temperature
        S =
                new float[nb_strat_temporelles][Population.X][Population.Y][Dataset_EVOL.prof_potentielles.length]; // Salinite
//        Th =new float[nb_strat_temporelles][Population.X][Population.Y]; // prof de melange / thermocline

        if (Ponte_Chla_SeaWiFS) {
            Chla = new double[12][Population.X][Population.Y];
        }


        //        strates = new int[Population.Pop_rec.size()][nb_strat_temporelles][X][Y][Dataset_EVOL.prof_potentielles.length];

        double[] TS; // vecteur de T, S, Th
        double z_grid; // pour la conversion en coordonnées de grille
        double ii, jj;
        int jour, d, i, j, k;
//        Nb_strates_date = new int[nb_strat_temporelles];
//        Nb_strates_poisson = new int[Population.Pop_rec.size()];



//        for (int p = 0; p < Population.Pop_rec.size(); p++) {
//            Nb_strates_poisson[p] = 0;
//        }


// ICI AMELIORATION POSSIBLE DE LA VITESSE EN UTILISANT "Java.util.map" plutôt que des tableau...
        // A voir, c'est pas sur...
        int month = 0;
        int new_month = 0;

        for (d = 0; d
                < nb_strat_temporelles; d++) {
            jour = d * resolution_temporelle_scanEnv;

            if (Ponte_Chla_SeaWiFS) {
                //System.out.println("month = " + month);

                new_month = (1 + (jour / 30) % 12);
                //System.out.println("new_month = " + new_month);

                if (month < new_month) {

                    //System.out.println("lecture CHLA du mois " + month);
                    try {
                        GetChla_SeaWiFSclim.load_SeaWiFS(jour);
                    } catch (IOException e) {
                        System.out.println("probleme dans ponte_init : " + e.getMessage());
                    }
                }
            }

            //System.out.println("Scan de l'environnement : jour  : " + jour);
            try {
                Dataset_EVOL.load_data(jour, year);
            } catch (IOException e) {
                System.out.println("probleme dans scan_env Chla 1 : " + e.getMessage());
            }

            for (i = 0; i
                    < Population.X; i++) {
                for (j = 0; j
                        < Population.Y; j++) {
                    ii = i + Population.x_min;
                    jj = j + Population.y_min;
                    if (Dataset_EVOL.isInWater(ii, jj)) {

                        if (flagGREG_dens) {
                            population.D_strates_spatiotemp[d][i][j] = 0;
                            population.D_strates_spatiotemp_finale[d][i][j] = 0;
                        }

// Lecture de la première profondeur k=0: on note T, S et Th (prof thermocline)
                        z_grid = Dataset_EVOL.depth2z(ii, jj, Dataset_EVOL.prof_potentielles[0]);
                        try {
                            TS = Dataset_EVOL.getFields_SaltTemp(ii, jj, z_grid, 0);
                            T[d][i][j][0] = (float) TS[0];// temp_init = evt[3];
                            S[d][i][j][0] = (float) TS[1];// salt_init =evt[4];
//                            Th[d][i][j] = (float) TS[2];

                        } catch (IOException e) {
                            System.out.println("probleme dans scan_env : " + e.getMessage());
                        }

                        if (Ponte_Chla_SeaWiFS) {
// Si c'est un nouveau mois, lecture de la Chla :
                            if (month < new_month) {
                                try {
                                    Chla[month][i][j] = GetChla_SeaWiFSclim.getCHLASeaWiFS(ii, jj);
                                } catch (ArrayIndexOutOfBoundsException e) {
                                    System.out.println("probleme dans scan_env Chla 2: " + e.getMessage());
                                }
                            }
                        }

// Lecture des autres prof : on note juste T et S                    
                        for (k = 1; k
                                < Dataset_EVOL.prof_potentielles.length; k++) {
                            z_grid = Dataset_EVOL.depth2z(ii, jj, Dataset_EVOL.prof_potentielles[k]);
                            try {
                                TS = Dataset_EVOL.getFields_SaltTemp(ii, jj, z_grid, 0);
                                T[d][i][j][k] = (float) TS[0];// temp_init = evt[3];
                                S[d][i][j][k] = (float) TS[1];// salt_init =evt[4];
                            } catch (IOException e) {
                                System.out.println("probleme dans scan_env : " + e.getMessage());
                            }

                        }
                    }
                }
            }
            month = 1 + (jour / 30) % 12; // ON CONSIDERE QUE LES MOIS ONT 30 JOURS...
        }
    }
//------------------------------------------------------------------------------ :   @

//---------------------------- ECRITURE DES RESULTATS (ASCII) ------------------ :   @
    public void closeResults() {

        System.out.println("on ferme le fihier ascii de sortie");
        try {
            resultsWriter.close();
        } catch (Exception e) {
            System.out.println("probleme en fermeture");
            System.exit(0);
        }
    }

// -----------------------------------------------------------------------------
    public void openResults(String strategy) {
        String filename;
            filename = output_dir + strategy + "_Ponte_Y" + compteur_year + ".txt";


        System.out.println(" Simulatio.829 : Création du fichier de sortie " + filename);

        try {
            resultsWriter = new FileWriter(filename);
        } catch (Exception e) {
            System.out.println("probleme en ouverture");
            System.exit(0);
        }

        boolean nom_des_colonnes = false;
        if (nom_des_colonnes) {
            //  String seoln = new String(new char[]{'\n'});
            StringBuffer s = new StringBuffer("year");
            s.append(";Month");
            s.append(";lat");
            s.append(";lon");
            s.append(";depth");
            s.append(";temp");
            s.append(";salt");
            s.append(";profthermo");
            s.append(";bathy");
            s.append(";phyto");
            s.append(";recruited?");
            s.append(";lat_final");
            s.append(";lon_final");
            s.append(";depth_final");
            //bathy finale
            //egg density

            s.append("\n");

            String ss = s.toString();
            try {
                resultsWriter.write(ss, 0, ss.length());
            } catch (Exception e) {
                System.out.println("probleme en ecriture");
                System.exit(0);
            }
        }
    }


//------------------------------------------------------------------------------
    public void writeResults() {

// Si pop est plein, on anregistre les donnees de Pop (ponte)
// Si pop est vide, alors on enregistre les donnees de Pop_rec (Rescrutement)	  

        String seoln = new String(new char[]{'\n'});
        for (int p = 0; p
                < population.Pop.size(); p++) {
            Poissons Po = (Poissons) population.Pop.get(p);
//    int mois = Po.day_init/15 +1;

            StringBuffer s = new StringBuffer(" " + year);
//        ss.append(";" + mois);
            s.append(";" + Po.day_init);
            s.append(";" + arondi(Po.lat_init)); // #3
            s.append(";" + arondi(Po.lon_init)); // #4
            s.append(";" + arondi(Po.depth_init));//#5
            s.append(";" + arondi(Po.temperature_natal));
            s.append(";" + arondi(Po.salinite_natal));
            s.append(";" + arondi(Po.profThermo_natal));
            s.append(";" + (int) Po.bathy_init); //#9
            s.append(";" + arondi(Po.phyto_init));
            int rec = Po.isRecruited ? 1 : 0;
            s.append(";" + rec); //#11
            s.append(";" + arondi(Po.lat));//#12
            s.append(";" + arondi(Po.lon));//#13
            s.append(";" + arondi(Po.depth));//#14
            s.append(";" + arondi(Po.poids));//#15 // densite de ponte obs correspondant a ce lieu de ponte (Si PONTE_OBS = true, sinon -999)
            s.append(";" + Po.egg_density);
            s.append(";" + Po.Body_length);//#17
            int retention = Po.isRetained_on_shelf ? 1 : 0;
            s.append(";" + retention);//#18
            //int nondisp = Po.isRecruitedGREGdens ? 1 : 0; // avant le 16 nov 2009: isRecruitedGREG
            //ss.append(";" + nondisp);//#19
            s.append(";" + Po.bathyActuelle);//#19 // remplace "nondis" (10/06/2010)
            // Tolerences :
            s.append(";" + Po.rayon_exploration_spatial);//#20
            s.append(";" + Po.rayon_exploration_temporel);//#21
            if (Strategie.equals("HE")) {
                s.append(";" + Po.tolerance_temperature);//#22
                s.append(";" + Po.tolerance_salinite);//#23
                s.append(";" + Po.tolerance_batymetrie);//#24
                s.append(";" + Po.tolerance_HBL);//#25
            }
            if (Strategie.equals("HS")) {
                s.append(";" + Po.lyaponov);//#22
            }
            if (Strategie.equals("BH")) {
                s.append(";" + Po.tempmax);//#22 // C'est la temperature maximale rencontrée
                s.append(";" + Po.tempmin);//#23 // C'est la temperature maximale rencontrée
            }

            s.append(seoln);
            String ss = s.toString();

            try {
                resultsWriter.write(ss, 0, ss.length());
            } catch (Exception e) {
                System.out.println("probleme en ecriture");
                System.exit(0);
            }
        }
    }
//-----------------------------------------------------------------------------

    public void openParam() {
        try {
            paramWriter = new FileWriter(output_dir + "Param.txt");
        } catch (Exception e) {
            System.out.println("probleme en ouverture" + output_dir);
            System.exit(0);
        }

        StringBuffer s = new StringBuffer("Parametres fixes pour le recrutement et les differentes strategies de ponte: \n");
        s.append("\n");

        s.append("A - PARAMETRE DE RECRUTEMENT : \n");
        s.append("ageAtDeath : " + (ageMinAtRecruitment + 1) + "\n");
        s.append("ageMinAtRecruitment" + ageMinAtRecruitment + "\n");
//	  ss.append("stayInRecDuration : " +stayInRecDuration +"\n");
        s.append("temp_lethal_min " + temp_lethal_min + " ; temp_lethal_max : " + temp_lethal_max + "\n");
        s.append("flagGREG_dens : " + flagGREG_dens + "\n");
        s.append("seuil_concentration : " + seuil_concentration + "\n");
        s.append("flagRET : " + flagRET_shelf + "\n");
        s.append("prof talu : " + Dataset_EVOL.prof_talu + "\n");
        s.append("flagRET_Chla : " + flagRET_Chla_SeaWiFS + "\n");                    // Recrutement si position dans la zone ou Chla SeaWiFS > Seuil_rec_Chla_SeaWiFS
        s.append("Seuil_rec_Chla_SeaWiFS : " + Seuil_rec_Chla_SeaWiFS + "\n");
        s.append("flagRET_CarbZ_PISCES : " + flagRET_CarbZ_PISCES + "\n");                    // Recrutement si position dans la zone ou Chla SeaWiFS > Seuil_rec_Chla_SeaWiFS
        s.append("Seuil_rec_CarbZ_PISCES : " + Seuil_rec_CarbZ_PISCES + "\n");


        s.append("B - PARAMETRE DE PONTE INITIALE : \n");
        s.append("nbFishesMois : " + nbFishesJour + "\n");
        s.append("prof_ponte_min : " + Dataset_EVOL.prof_ponte_min + " ; prof_ponte_max : " + Dataset_EVOL.prof_ponte_max + "\n");
        s.append("lat_min, lat_max, lon_min, lon_max : " + Dataset_EVOL.lat_min + " , " + Dataset_EVOL.lat_max + " , " + Dataset_EVOL.lon_min + " , " + Dataset_EVOL.lon_max + "\n");
        s.append("bathy_max : " + Dataset_EVOL.bathy_max + "\n");

        s.append("\n");

        s.append("C - PARAMETRE DE STRATEGIE DE PONTE : \n");
        s.append("patchiness_rad_HS (en Km): " + patchiness_rad_HS + "\n");
        s.append("time_rad_HS : " + time_rad_HS + "\n");
//        ss.append("patchiness_rad_HE (en Km) : " + patchiness_rad_HE + "\n");
//        ss.append("rayon_exploHE (en Km): " + rayon_exploHE + "\n");
        s.append("patchiness_rad_OPP (en Km): " + patchiness_rad_OPP + "\n");
        s.append("Ponte_Chla_SeaWiFS : " + Ponte_Chla_SeaWiFS + "\n");
        s.append("Seuil_ponte_Chla_SeaWiFS = " + Seuil_ponte_Chla_SeaWiFS + "\n");

        s.append("Range_rayon_exploration_spatial max : " + rayon_exploration_spatial_max + "\n");
        s.append("Range_rayon_exploration_spatial min : " + rayon_exploration_spatial_min + "\n");
        s.append("Range_rayon_exploration_temporel max: " + rayon_exploration_temporel_max + "\n");
        s.append("Range_rayon_exploration_temporel min: " + rayon_exploration_temporel_min + "\n");
        s.append("Range_tolerance_temperature : " + Range_tolerance_temperature + "\n");
        s.append("Range_tolerance_salinite : " + Range_tolerance_salinite + "\n");
        s.append("Range_tolerance_batymetrie : " + Range_tolerance_batymetrie + "\n");
        s.append("Range_tolerance_HBL : " + Range_tolerance_HBL + "\n");
        s.append("flag_Densite_max : " + flag_Densite_max + "\n");
        s.append("D_max : " + Population.D_max_spatiotemp + "\n");


        s.append("HomingGeographique = " + HomingGeographique + "\n");
        s.append("HomingTemporel = " + HomingTemporel + "\n");
        s.append("HomingTemperature = " + HomingTemperature + "\n");
        s.append("HomingSalinite = " + HomingSalinite + "\n");
        s.append("HomingBathymetrie = " + HomingBathymetrie + "\n");
        s.append("HomingHBL = " + HomingHBL + "\n");
        s.append("\n");

        /*
        ss.append("temp_marge : " + temp_marge + "\n");
        ss.append("salt_marge : " + salt_marge + "\n");
        ss.append("h_marge : " + h_prop + "\n");
        ss.append("Th_marge : " + Th_marge + "\n");
        ss.append("prof_marge : " + prof_marge + "\n");
         */
        s.append("E - PARAMETRES INDIVIDUS : \n");

        s.append("egg_buoyancy_flag = " + flagBuoy + "\n"); //[boolean] flotabilité des oeuf: si false neutre, si true
        s.append("egg_buoyancy = " + egg_excess_buoyancy + "\n"); // d'après Gorant et al. 2007
        s.append("egg_duration = " + egg_duration + "\n");
        s.append("flag_Growth = " + flagGrowth + "\n"); //[boolean]
        s.append("flag_Swim_O2 = " + flagSwim_O2 + "\n"); //[boolean]
        s.append("flag_Swim_random = " + flagSwim_random + "\n"); //[boolean]
        s.append("Swim_upperlimit = " + Swim_upperlimit + "\n");
        s.append("Swim_lowerlimit = " + Swim_lowerlimit + "\n");
        s.append("\n");

        s.append("F - PARAMETRES DIVERS : \n");
        s.append("dt_advec" + dt_advec/60 + "heures \n");
        s.append("Topp, Sopp :" + Topp + " ; " + Sopp + "\n");
        s.append("Prof de thermocline Opp :" + Th_opp + "\n");
        s.append("nbre_generations : " + nbre_generations + "\n");
        s.append("nb_jours_par_ans : " + nb_jours_par_ans + "\n");

        s.append("\n");
        String ss = s.toString();

        try {
            paramWriter.write(ss, 0, ss.length());
        } catch (Exception e) {
            System.out.println("probleme en ecriture");
            System.exit(0);
        }
    }
//-----------------------------------------------------------------------------

    public void closeParam() {
        System.out.println("on ferme");
        try {
            paramWriter.close();
        } catch (Exception e) {
            System.out.println("probleme en fermeture");
            System.exit(0);
        }
    }
// -----------------------------------------------------------------------------

    public float arondi(double nombre_long) {
        float nombre_cour = (float) Math.round(nombre_long * 100) / 100;
        return (nombre_cour);
    }
// -----------------------------------------------------------------------------


    void TEST_DEB(){
        // 1er lecture des champs dynamique
        try {
            dataset_EVOL.init();
        } catch (IOException e) {
            System.out.println("blem 1: " + e.getMessage());
            System.exit(0);
        }
population.ponte_init_TEST_DEB();

        // Calcul du nombre d'iteration a faire lors de l'advections des individus (OEUFS ET LARVES):
        nbiter = (int) (((double) Dataset_EVOL.dt_HyMo * 24*60) / (double) dt_advec); //(dt_HyMo est en jour et dt_advec en heures)
        dt = 1.0f / nbiter;
        System.out.println("Pas de temps pour l'advection des larves : " + dt_advec + "heures");
        System.out.println("Soit " + nbiter + " interpolations entre les sorties roms.");
        // nbre d'iteration par pas de temps des sorties roms         // exemple : peru : sorties a 2 jours--> 12 iteration = une actualisation de champs de vitesse toute les 4 heure

                while (compteur_year < nbre_generations) {
            year = Dataset_EVOL.yearlist100[compteur_year];

            System.out.println("<<<<<<<<<<<<<<<< TEST_DEB - GENERATION " + compteur_year + ">>>>>>>>>>>>>>>>");

        // CREATION DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" CREATION DU NetCDF TEST  " );
            try {
            netcdf_outpout_writer.openResults_NetCDF();
            netcdf_outpout_writer.open_Biom2D_NetCDF();
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf TEST_DEB : " + ex.getMessage());
            }
}
        compteur_t_rec_intra_year =0;
                population.stepSerial_TEST_DEB();

                // FERMETURE DU NetCDF --------------------------------------------------
if(outpout_netcdf_flag){
        System.out.println(" FERMETURE DE L'ANCIEN  NetCDF TEST_DEB " );
            try {
            netcdf_outpout_writer.CloseResults_NetCDF();
            netcdf_outpout_writer.Close_Biom2D_NetCDF();
            } catch (Exception ex) {
                System.out.println("probleme dans creation netcdf TEST_DEB : " + ex.getMessage());
            }
}

            compteur_year = compteur_year + 1;
        } // Fin du while (compteur_year < nbre_generations)
        if (compteur_year == nbre_generations) { // en principe ça devrait être le cas!!
            System.out.println(" c'est fini pour : " + Strategie);
            System.exit(0);
        }
    }
}
