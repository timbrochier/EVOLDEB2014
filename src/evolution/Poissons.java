package evolution;

//import java.awt.*;
//import java.util.*;
import java.io.*;
import static java.lang.Double.*;

//import evolution.DebParticleLayer;
//import evolution.Kinesis;
import java.util.Random;

public class Poissons {
// ----------------------------------------------------------------------------

    int compteur_advec = 0;
//****************************************
    public double x, y, z, temp, salt, phyto, lat, lon,
            lyaponov, grad_temp;
    
        public float Body_length, depth, bathyActuelle;
    public long id;
    public int age;
    int compteur_phyto, ss, Nb_strates;
    public boolean living, isRecruited, isRetained_on_shelf, isRecruitedGREG,
            isRecruitedGREGdens, isdead_predation, isdead_temp, first_spawn;
    // Variables de mémoire : (enrgistrée une fois au début pour chaque poisson)
    // 1 - environnement natal (ou initial)
    public double lat_init, lon_init, depth_init, profThermo_natal, DistVoisin_init,
            DistVoisin_finale, x_init, y_init, z_init, temperature_natal, salinite_natal,
            bathy_init, phyto_init, egg_density;
    // Variables de test recrutement : (enrgistrée une fois à la fin pour chaque poisson)
    double chla_fin, Boufe_avant;
    // 3 - tolerances environnementales, spatiales et temporelles (hereditaire)
    double tolerance_temperature, tolerance_salinite, tolerance_HBL;
    public int day_init, year_init; 
    int nbre_voisins_ini, nbre_voisins, rayon_exploration_spatial,
            rayon_exploration_temporel, tolerance_batymetrie;
    double[] uv, TS, plankton;
    int[] Tolerances_spatiotemps;
    double[] Tolerances_TS;
    boolean[][][] strates;
    // Variables de super individus :
    double poids;
    // Pour enregistrer les temperatures max et min rencontree par chaque
    // individu (6 octobre 2011): 
    float tempmax, tempmin;
    public double bouffe;
    double bouffe_max_0, bouffe_max_1, bouffe_max_2, bouffe_max_3;
    // bouffe_max_0 : memoir de la plus forte concentration de plankton[0] vue
    //ArrayList positions_OK;
//    private final List<Position_ponte>
    public DebParticleLayer DEB;
    public double Q, Qold; // Habitat quality
    float CarCapa; // Habitat Carrying Capacity
    float densite_biomasse_locale_par_km2_autres_super_individu, densite_biomasse_de_ce_banc_par_km2, densite_biomasse_locale_par_km2_TOTALE;

    ; // Densite dans la zone
    public double S; // Nombre d'individu que représente ce super-individu
    double weight, weight_old; // Poids en grammes
    public double[] V_kinesis, V_kinesis_old; // composante x,y,z de la nage donne par la kinesis en m/s
    double dx_kinesis_swim_km, dy_kinesis_swim_km, dz_swim_m;// composante x,y,z de la nage donne par la kinesis en km
    double f_reponse_fonctionnel;
    static float dist_max_par_dt, dist_parcourue_km;
    int compteur_standstill; // pour enregistrer le nombre de fois d'affiler que l'indiv est collé à la côte
    // --> si >  100 fois on considere l'indiv comme perdu (échoué)
    public long Nb_eggs_spawned_dt, Nb_SI_spawned_dt;

// ----------------------------------------------------------------------------
    public Poissons(int jour) {
        ss = 0;
        boolean outZones = true;
        while (outZones) {

            int bat_max = Dataset_EVOL.bathy_max;

            x = Population.x_min + Math.random() * (Population.x_max - 1 - Population.x_min);
            y = Population.y_min + Math.random() * (Population.y_max - 1 - Population.y_min);
            poids = -999;

            // if ((x < Population.x_min) || (y < Population.y_min)) {
            //     System.out.println(" x = " + x + " y = " + y + " depth = " + depth);
            //}

            int ix = (int) x;
            int iy = (int) y;
            double xi = x;
            double yi = y;
            depth = (float) -(Dataset_EVOL.prof_ponte_min + Math.random() * (Dataset_EVOL.prof_ponte_max - Dataset_EVOL.prof_ponte_min));

            z = Dataset_EVOL.depth2z(x, y, depth - Simulation.decalage_zeta);

            outZones = (!Dataset_EVOL.isInWater(ix, iy))
                    || (depth < getDepthBottom(x, y))
                    || (x > (Dataset_EVOL.nx - 2.0f))
                    || (x < 2.0f)
                    || (y > (Dataset_EVOL.ny - 3.0f))
                    || (y < 2.0f)
                    || (Math.abs(x - xi) > 1.0f)
                    || (Math.abs(y - yi) > 1.0f)
                    || (Dataset_EVOL.isCloseToCost(x, y))
                    || (Dataset_EVOL.getBathy(ix, iy) > bat_max);


            ss++;
            if (ss > 10000) {
                System.out.println(" impossible  (1) " + Population.Pop.size());
                System.out.println(" depth = " + depth + " ;x= " + x + " ;y= " + y);
                System.out.println("(!Dataset_EVOL.isInWater(x, y))  = " + (!Dataset_EVOL.isInWater(x, y)));
                System.out.println("(x > (Dataset_EVOL.nx - 3.0f)) = " + (x > (Dataset_EVOL.nx - 3.0f)));
                System.out.println("(Dataset_EVOL.isCloseToCost(x, y)) = " + (Dataset_EVOL.isCloseToCost(x, y)));
                System.out.println("(depth < getDepthBottom(x, y)) = " + (depth < getDepthBottom(x, y)));
                System.out.println("depth = " + depth);
                System.out.println("getDepthBottom(x, y) = " + getDepthBottom(x, y));
                System.out.println("(Dataset_EVOL.getBathy(ix, iy) > Dataset_EVOL.bathy_max) = " + (Dataset_EVOL.getBathy(ix, iy) > Dataset_EVOL.bathy_max));
                System.out.println("(Dataset_EVOL.getBathy(ix, iy) = " + Dataset_EVOL.getBathy(ix, iy));
                System.out.println("(Dataset_EVOL.bathy_max = " + Dataset_EVOL.bathy_max);

                System.out.println(" x = " + x + " y = " + y + " depth = " + depth);
                System.exit(0);
            }
        }
        if (Simulation.flag_Densite_max || Simulation.flagGREG_dens) {
            int ii = Math.round(Math.round(x));
            int jj = Math.round(Math.round(y));
            int i = ii - Population.i_min;
            int j = jj - Population.j_min;
            int d = Math.round(jour / Simulation.resolution_temporelle_scanEnv);
            Population.D_strates_spatiotemp[d][i][j]++;
//System.out.println("222 jour = " + jour + " ; d= " + d + "  ; Simulation.resolution_temporelle_scanEnv = " + Simulation.resolution_temporelle_scanEnv);
        }
        init_tolerances();
        init(jour);
        get_identity();
    }

    // ----------------------------------------------------------------------------
    public Poissons(double fx, double fy, double fdepth, int jour, double radius, int[] Heritage_tolerances_spatiotemps, double[] Heritage_tolerances_TS) {

        ss = 0;
        double p = Simulation.vertical_patchiness_rad; // rayon sur la profondeur = +-2m     
        boolean outZones = true;
        boolean Densite_trop_importante = false;
        while (outZones) {
            double r, t;
            x = -1.f;
            y = -1.f;
            while ((x < Population.x_min) || // ajout 30 septembre; 5 Nov: modif "Population.sponge_grid"
                    (x > Population.x_max)
                    || (y < Population.y_min) || // ajout 30 septembre
                    (y > Population.y_max)) {
                t = Math.random() * 2.0f * Math.PI; // direction du shift
                r = Math.random() * radius;         // radius  = rayon de patchiness (en nbre de mailles)                
                x = fx + r * Math.cos(t);	        // distance horiz x
                y = fy + r * Math.sin(t);	        // distance horiz y

                //System.out.println("x_min = " + Population.x_min +" ; x_max = " + Population.x_max + "  -> x = " + x);
                //System.out.println("y_min = " + Population.y_min +" ; y_max = " + Population.y_max + "  -> y = " + y);

            }
            // ajout 30 septembre

            //int ix = (int) x;
            //int iy = (int) y;
            depth = (float) (fdepth + 2 * p * (Math.random() - 0.5f)); // shift sur la profondeur de ponde
            // Securite por pas qu'elles pondent en l'air :
            bathyActuelle = (float) getDepthBottom(x, y);

            if (depth > -0.5) {
                depth = -1;
            }
            if (depth < bathyActuelle) {
                depth = (bathyActuelle + Math.abs(bathyActuelle / 2));
            }
            int ix = (int) x;
            int iy = (int) y;

            outZones = ((!Dataset_EVOL.isInWater(ix, iy))
                    || (x > (Dataset_EVOL.nx - 2.0f))
                    || (x < 2.0f)
                    || (y > (Dataset_EVOL.ny - 3.0f))
                    || (y < 2.0f))
                    || (x < Population.x_min) || // ajout 30 septembre
                    (x > Population.x_max) || // ajout 30 septembre
                    (y < Population.y_min) || // ajout 30 septembre
                    (y > Population.y_max) || // ajout 30 septembre
                    (Dataset_EVOL.isCloseToCost(x, y))
                    || (depth < getDepthBottom(x, y));
            ss++;
            if (ss > 7000) {
                System.out.println(" impossible (2) ");
                System.out.println(" depth = " + depth + " ;x= " + x + " ;y= " + y);
                System.out.println("(!Dataset_EVOL.isInWater(x, y))  = " + (!Dataset_EVOL.isInWater(x, y)));
                System.out.println("(x > (Dataset_EVOL.nx - 3.0f)) = " + (x > (Dataset_EVOL.nx - 3.0f)));
                System.out.println("(Dataset_EVOL.isCloseToCost(x, y)) = " + (Dataset_EVOL.isCloseToCost(x, y)));
                System.out.println("(depth < getDepthBottom(x, y)) = " + (depth < getDepthBottom(x, y)));
                System.out.println("depth = " + depth);
                System.out.println("getDepthBottom(x, y) = " + getDepthBottom(x, y));
                System.out.println(Math.round(jour / Simulation.resolution_temporelle_scanEnv));
                System.out.println("Densite_trop_importante : " + Densite_trop_importante);
                x = fx;
                y = fy;
                z = Dataset_EVOL.depth2z(x, y, fdepth - Simulation.decalage_zeta);
                System.exit(0);
                outZones = false;
            }
        }
        z = Dataset_EVOL.depth2z(x, y, depth - Simulation.decalage_zeta);

// ON MEMORISE L'ENVIRONNEMENT DE NAISSANCE :
        day_init = jour;

// HERITAGE DES TOLERANCES :
// Memoire des tolerances des parents, avec une mutation aleatoire :
        double alea;
        int i = 2;
        double R_alea;

        // Rayon d'exploration spatial :----------------------------------------
        if (Simulation.HomingGeographique) {
            if (Simulation.rayon_exploration_spatial_min == Simulation.rayon_exploration_spatial_max) {
                rayon_exploration_spatial = (int) Simulation.rayon_exploration_spatial_min;
            } else {
                R_alea = Population.mean_grid_size;
                alea = (2 * Math.random() * R_alea) - R_alea;
                if (Heritage_tolerances_spatiotemps[0] + alea > Simulation.rayon_exploration_spatial_max) {
                    rayon_exploration_spatial = (int) (Simulation.rayon_exploration_spatial_max * 2 - Heritage_tolerances_spatiotemps[0] - ((int) Math.round(alea)));
                } else {
                    rayon_exploration_spatial = (int) Simulation.rayon_exploration_spatial_min + ((int) Math.abs(Heritage_tolerances_spatiotemps[0] - Simulation.rayon_exploration_spatial_min + alea));
                }
            }
        } else {
            rayon_exploration_spatial = 9999999;
        }

        //----------------------------------------------------------------------
        // Rayon d'exploration temporel : 
        if (Simulation.HomingTemporel) {
            if (Simulation.rayon_exploration_temporel_min == Simulation.rayon_exploration_temporel_max) {
                rayon_exploration_temporel = (int) Simulation.rayon_exploration_temporel_min;
            } else {
                R_alea = Math.max(Simulation.resolution_temporelle_scanEnv, Dataset_EVOL.dt_HyMo);
                alea = 2 * Math.random() * R_alea - R_alea;
                if (Heritage_tolerances_spatiotemps[1] + alea > Simulation.rayon_exploration_temporel_max) {
                    rayon_exploration_temporel = (int) Math.max(Simulation.rayon_exploration_temporel_max * 2 - Heritage_tolerances_spatiotemps[1] - ((int) Math.round(alea)), Simulation.rayon_exploration_temporel_min);
                } else {
                    rayon_exploration_temporel = (int) Simulation.rayon_exploration_temporel_min + (int) Math.abs(Heritage_tolerances_spatiotemps[1] - Simulation.rayon_exploration_temporel_min + alea);
                }
            }
        } else {
            rayon_exploration_temporel = 183;
        }
        //----------------------------------------------------------------------
        // Tolérance aux variations de Bathymetrie :
        if (Simulation.HomingBathymetrie) {
            R_alea = 10;
            alea = 2 * Math.random() * R_alea - R_alea;
            if (Heritage_tolerances_spatiotemps[2] + alea > Simulation.Range_tolerance_batymetrie) {
                tolerance_batymetrie = (int) Math.round((Simulation.Range_tolerance_batymetrie * 2 - Heritage_tolerances_spatiotemps[2] - alea));
            } else {
                tolerance_batymetrie = (int) Math.round(Math.abs(Heritage_tolerances_spatiotemps[2] + alea));
            }
        } else {
            tolerance_batymetrie = 99999;
        }
        //----------------------------------------------------------------------
        if (Simulation.HomingTemperature) {
            // Tolérance aux variation de température : mutation de + ou - 0.5
            alea = (float) (Math.random() - 0.5f) / 10;
            if (Heritage_tolerances_TS[0] + alea > Simulation.Range_tolerance_temperature) {
                tolerance_temperature = Simulation.Range_tolerance_temperature * 2.0f - Heritage_tolerances_TS[0] - alea;
            } else {
                tolerance_temperature = Math.abs(Heritage_tolerances_TS[0] + alea);
            }
        } else {
            tolerance_temperature = 99999;
        }
        //----------------------------------------------------------------------
        if (Simulation.HomingSalinite) {
            // Tolérance aux variation de salinité : mutation de + ou - 0.1 PSU
            alea = (float) (0.02f * Math.random() - 0.01);
            if (Heritage_tolerances_TS[1] + alea > Simulation.Range_tolerance_salinite) {
                tolerance_salinite = Simulation.Range_tolerance_salinite * 2.0f - Heritage_tolerances_TS[1] - alea;
            } else {
                tolerance_salinite = Math.abs(Heritage_tolerances_TS[1] + alea);
            }
        } else {
            tolerance_salinite = 99999;
        }

        //----------------------------------------------------------------------
        // Tolérance aux variation de HBL : mutation de + ou - 5m
        if (Simulation.HomingHBL) {
            alea = (float) (5.0f * Math.random() - 2.5f);
            if (Heritage_tolerances_TS[2] + alea > Simulation.Range_tolerance_HBL) {
                tolerance_HBL = Simulation.Range_tolerance_HBL * 2.0f - Heritage_tolerances_TS[2] - alea;
            } else {
                tolerance_HBL = Math.abs(Heritage_tolerances_TS[2] + alea);
            }
        } else {
            tolerance_HBL = 99999;
        }

        // On remplit les tableaux "Tolerances_spatiotemps" et "Tolerances_TS"
        // avec les nouvelles valeurs "mutées"
        Tolerances_spatiotemps = new int[]{rayon_exploration_spatial, rayon_exploration_temporel, tolerance_batymetrie};
        Tolerances_TS = new double[]{tolerance_temperature, tolerance_salinite, tolerance_HBL};
//        id = Population.compteur_id_poissons;
        get_identity();

    }

    // ----------------------------------------------------------------------------
    public Poissons(double fx, double fy, double fdepth, int jour, double radius, int Nbeggs) {

        S = Nbeggs;

        ss = 0; // compteur pour la boucle while ci-dessous
        double p = Simulation.vertical_patchiness_rad; // rayon sur la profondeur = +-2m
        double r, t;
        boolean outZones = true;
        boolean Densite_trop_importante = false;
        while (outZones) {
            //x = -1.f;
            //y = -1.f;
            t = Math.random() * 2.0f * Math.PI; // direction du shift
            r = Math.random() * radius;         // radius  = rayon de patchiness (en nbre de mailles)
            x = fx + r * Math.cos(t);	        // distance horiz x
            y = fy + r * Math.sin(t);	        // distance horiz y

            /*
            if ((x < Population.x_min) || // ajout 30 septembre; 5 Nov: modif "Population.sponge_grid"
            (x > Population.x_max)
            || (y < Population.y_min) || // ajout 30 septembre
            (y > Population.y_max)) {

            System.out.println(" PONTE HORS DOMAINE");
            System.out.println("x_min = " + Population.x_min + " ; x_max = " + Population.x_max + "  -> x = " + x);
            System.out.println("y_min = " + Population.y_min + " ; y_max = " + Population.y_max + "  -> y = " + y);
            }
             *
             */
            // ajout 30 septembre

            //int ix = (int) x;
            //int iy = (int) y;
            depth = (float) (fdepth + 2 * p * (Math.random() - 0.5f)); // shift sur la profondeur de ponde
            // Securite por pas qu'elles pondent en l'air :
            bathyActuelle = (float) getDepthBottom(x, y);
            if (depth > -0.5) {
                depth = -1;
            }
            if (depth < bathyActuelle) {
                depth = (float) (bathyActuelle + Math.abs(bathyActuelle / 2));
            }
            int ix = (int) x;
            int iy = (int) y;

            outZones = ((!Dataset_EVOL.isInWater(ix, iy))
                    || (x > (Dataset_EVOL.nx - 2.0f))
                    || (x < 2.0f)
                    || (y > (Dataset_EVOL.ny - 3.0f))
                    || (y < 2.0f))
                    || (x < Population.x_min) || // ajout 30 septembre
                    (x > Population.x_max) || // ajout 30 septembre
                    (y < Population.y_min) || // ajout 30 septembre
                    (y > Population.y_max) || // ajout 30 septembre
                    (Dataset_EVOL.isCloseToCost(x, y))
                    || (depth < getDepthBottom(x, y));
            ss++;
            if (ss > 7000) {
                System.out.println(" impossible (3) ");
                System.out.println(" depth = " + depth + " ;x= " + x + " ;y= " + y);
                System.out.println("(!Dataset_EVOL.isInWater(x, y))  = " + (!Dataset_EVOL.isInWater(x, y)));
                System.out.println("(x > (Dataset_EVOL.nx - 3.0f)) = " + (x > (Dataset_EVOL.nx - 3.0f)));
                System.out.println("(Dataset_EVOL.isCloseToCost(x, y)) = " + (Dataset_EVOL.isCloseToCost(x, y)));
                System.out.println("(depth < getDepthBottom(x, y)) = " + (depth < getDepthBottom(x, y)));
                System.out.println("depth = " + depth);
                System.out.println("getDepthBottom(x, y) = " + getDepthBottom(x, y));
                System.out.println(Math.round(jour / Simulation.resolution_temporelle_scanEnv));
                System.out.println("Densite_trop_importante : " + Densite_trop_importante);
                System.out.println("Point de départ fx = " + fx + " , fy = " + fy + " fdepth = " + fdepth);
                System.exit(0);

                x = fx;
                y = fy;
                z = Dataset_EVOL.depth2z(x, y, fdepth - Simulation.decalage_zeta);
                outZones = false;
            }
        }
        z = Dataset_EVOL.depth2z(x, y, depth - Simulation.decalage_zeta);

// ON MEMORISE L'ENVIRONNEMENT DE NAISSANCE :
        init(jour);
        // = Population.compteur_id_poissons;
        get_identity();

        //        System.out.println("ID poisson pondu = " + id);
    }
//--------------------------------------------------------------------------

    void init_tolerances() {
        rayon_exploration_spatial = Math.round(Math.round(Simulation.rayon_exploration_spatial_max * Math.random()));
        rayon_exploration_temporel =
                Math.round(Math.round(Simulation.rayon_exploration_temporel_max * Math.random())) + 1;
        tolerance_batymetrie =
                Math.round(Math.round(Simulation.Range_tolerance_batymetrie * Math.random()));

        tolerance_temperature =
                Simulation.Range_tolerance_temperature * Math.random();
        tolerance_salinite =
                Simulation.Range_tolerance_salinite * Math.random();
        tolerance_HBL =
                Simulation.Range_tolerance_HBL * Math.random();

        // stockage dans "Tolerance" (pour n'avoir qu'une variable à passer aux heritiers)
        /*
        Tolerances[0] = rayon_exploration_spatial;
        Tolerances[1] = rayon_exploration_temporel;
        Tolerances[2] = tolerance_temperature;
        Tolerances[3] = tolerance_salinite;
        Tolerances[4] = tolerance_batymetrie;
        Tolerances[5] = tolerance_HBL;
         */
        Tolerances_spatiotemps =
                new int[]{rayon_exploration_spatial, rayon_exploration_temporel, tolerance_batymetrie};
        Tolerances_TS =
                new double[]{tolerance_temperature, tolerance_salinite, tolerance_batymetrie, tolerance_HBL};
    }

// -----------------------------------------------------------------------------
    void init(
            int jour) {

// DEB -- -- -- -- - - - --
        DEB = new DebParticleLayer();
        DEB.init();
        first_spawn = true;
// ---------------------------

        positGeog2D();  //(x,y) ==> (lat,lon)                
        day_init = jour;
        year_init = Simulation.year;
        living = true;
        isRecruited = false;
        isdead_temp = false;
        age = 0;
        Body_length = 0.1f; // cm
        weight = 0.01; //g

        compteur_standstill = 0;


        // (float) util.GrowthModel.LENGTH_INIT; // taille de l'oeuf

        // localisation initiale :
        lat_init = lat;
        lon_init = lon;
        depth_init = depth;
        bathy_init = getDepthBottom(x, y);
        dist_parcourue_km = 0;

        x_init = x;
        y_init = y;
        z_init = z;

        lyaponov = -999;

        // Pour enregistrement des temp max et min rencontrée (pour BH, 6 octobre 2011)
        tempmax = -999;
        tempmin = 999;

        DistVoisin_init = 999999999; // ( initialement mise a l'infinie, puis est calculee si on active le module greg)
        DistVoisin_finale = 999999999;
        compteur_phyto = 0;


        // Super indiv -------------
        S = (int) Simulation.nbeggs_max_per_SI; // Nbre dd'individus que représente le super indiv.
        // Calibrer ce nombre tq pour une mortalite naturelle on arrive a des bancs adulte de ~ 1 tonne
        //
        Q = 0; // init habitat quality
        V_kinesis_old = new double[]{0, 0, 0}; // composante x,y,z de la nage donne par la kinesis
        V_kinesis = new double[]{0, 0, 0}; // composante x,y,z de la nage donne par la kinesis

    }
    /*
    // -----------------------------------------------------------------------------
    void deplace_ichthyoplankton() {
    //System.out.println("NEW DEPLACE ... jj = " + jj + "****************************");
    int it_debut, it_fin;
    it_debut =
    (int) Population.nb_jours_ecoule_entre_roms_records * 24 / (int) Simulation.dt_advec;
    it_fin =
    (int) (Population.nb_jours_ecoule_entre_roms_records + 1) * 24 / (int) Simulation.dt_advec;
    // *****************************************************
    for (int it = it_debut; it < it_fin; it++) {

    double frac = it * Simulation.dt;
    try {
    // LECTURE/INTERPOLATION DES CHAMPS DE VITESSE U, V et W :
    uvw = Dataset_EVOL.getFields_uvw(x, y, z, frac);
    // LECTURE/INTERPOLATION DES CHAMPS TEMPERATURE ET SALINITE :
    TS = Dataset_EVOL.getFields_SaltTemp(x, y, z, frac);
    // LECTURE/INTERPOLATION DES CHAMPS PLANKTON ET OXYGENE :
    if (Simulation.flagGrowth) { //|| Simulation.flagSwim_O2) { // rétablir si lecture O2 Pisces necessaire pour swim (retablir dans Poissons et dans Dataset_EVOL)
    double[] pGrid = {x, y, z};
    //plankton = Dataset_EVOL.getPlankton_PISCES(pGrid, frac);
    plankton = Dataset_EVOL.getPlankton_PISCES(pGrid, frac);
    bouffe = plankton[1];

    CarCapa = CarCapa + CarCapa_locale(bouffe, Simulation.dt_advec);

    }

    // LECTURE/INTERPOLATION DES CHAMPS DE DEPLACEMENT DU A LA TURBULENCE SOUS-MAILLE
    if (Simulation.flagTurbulence) {
    // INTRODUIRE ICI LE MODULE TURBULENCE DE ICHTHYOP
    }

    temp = TS[0];
    salt = TS[1];

    // MODIF DU 28 juin 2013 : on ne fait pas le deplacement si celui-ci nous amene a terre :
    if (Dataset_EVOL.isInWater(x + uvw[0], y + uvw[1])) {
    x = x + uvw[0];
    y = y + uvw[1];
    z = z + uvw[2];
    }

    if ((Simulation.flagBuoy) || (Simulation.flagSwim_O2) || Simulation.flagSwim_random) {
    depth = (float) Dataset_EVOL.z2depth(x, y, z);
    depth = (float) (depth + Simulation.decalage_zeta); // Pour corriger le souscis avec le run preindus
    }

    bathyActuelle = getDepthBottom(x, y); //Dataset_EVOL.getDepth(x, y, 0);

    ///////////// ----- E X I T  O F  D O M A I N ----------------------------------
    auBord();
    if (!living) {
    break;
    }

    ///////////// ----- G R O W T H ------------------------------------------------
    if (Simulation.flagGrowth) {
    //double Diatoms = plankton[0];
    //double NanoPhyto = plankton[1];
    //double MicroZoo = plankton[2];
    //double MesoZoo = plankton[3];
    //length = util.GrowthModel.grow(length, temp, Diatoms, NanoPhyto, MicroZoo, MesoZoo);

    DEB.execute(temp, bouffe, (double) Simulation.dt_advec / 24);

    /*
    System.out.println("Diatoms = " + Diatoms);
    System.out.println("NanoPhyto = " + NanoPhyto);
    System.out.println("MicroZoo = " + MicroZoo);
    System.out.println("MesoZoo = " + MesoZoo);

    }

    ///////////// ----- S W I M ( D V M ) ------------------------------------------
    // MIGRATION VERTICALE (DVM SURFACE - Oxycline) :
    if (Simulation.flagSwim_random) {
    if (!Simulation.flagGrowth && (age > Simulation.egg_duration)) {
    //                            || Simulation.flagGrowth && (util.GrowthModel.getStage(length) == 2)) {
    // EN L'ABSENCE DE DONNEE PLUS PRECISE :
    depth = (float) (-Math.random() * Simulation.Swim_lowerlimit);
    }

    } else if (Simulation.flagSwim_O2) {
    if (!Simulation.flagGrowth && (age > Simulation.egg_duration)) {
    //                            || Simulation.flagGrowth && (util.GrowthModel.getStage(length) == 2)) {

    // En lisant la prof de l'oxycline 1ml.L-1 dans les données :
    double prof_Oxy = GetOxyclinDepth_BIDOUILLE.getOxyclin_depth(x, y);
    // Pour les cas ou la couche 1ml.L-1 n'est pas définie :
    //prof_Oxy = Math.max(Math.max(prof_Oxy, bathyActuelle), -100); // prof max de migration = 100m (si bathy et oxycline plus profonde)
    //System.out.println("Max = prof_Oxy = " + prof_Oxy + "  ; bathyActuelle" + bathyActuelle);
    depth = (float) (Math.random() * Math.max(prof_Oxy, -50)); // (prof_Oxy est négatif)
    }
    }
    } catch (IOException e) {
    System.out.println("youstone on a un probleme : " + e.getMessage());
    }

    // DEBUT Correction sur les valeurs aberrantes de depth : ----------------------
    ss = 0;
    boolean refaire_prof = false;
    while ((depth < bathyActuelle + 0.1) || (depth > -0.1)) {
    refaire_prof = true;
    ss++;

    if (depth > -0.1) {
    depth = (float) (-0.1 - Math.random() * 3);
    } else if (depth < bathyActuelle + 1) {
    depth = (float) (bathyActuelle + 2 + Math.random() * 3);
    }

    if (ss > 100) {
    System.out.println(" impossible  (3) - bathyActuelle = " + bathyActuelle);
    System.out.println(" x = " + x + " y = " + y + " depth = " + depth);
    System.exit(0);
    }
    }
    /// FIN de correction sur les valeurs aberrantes de depth ----------------------
    // Passage de la profondeur en coordonnées Grille ******************************
    if ((Simulation.flagBuoy) || (Simulation.flagSwim_O2) || refaire_prof) {
    depth = (float) (depth - Simulation.decalage_zeta); // Pour corriger le souscis avec le run preindusF
    z = Dataset_EVOL.depth2z(x, y, depth);
    }
    positGeog3D(); // Passage des coordonées grille aux coordonnée lon lat depth
    }
    }
     */
// -----------------------------------------------------------------------------

    void positGeog2D() {
        double[] po = Dataset_EVOL.grid2Geo(x, y);
        lat =
                po[0];
        lon =
                po[1];
        //depth = po[2];
    }
// -----------------------------------------------------------------------------

    void positGeog3D() {
        double[] po = Dataset_EVOL.grid2Geo(x, y, z);

        if ((lat == po[0]) && (lon == po[1])) {
            compteur_standstill = compteur_standstill + 1;
            //     System.out.println("compteur_standstill + 1 = " + compteur_standstill);
        }
        lat = po[0];
        lon = po[1];
        depth = (float) po[2];
    }
// -----------------------------------------------------------------------------
// Donne la profondeur (bathy) au point xRho, yRho.

    double getDepthBottom(double xRho, double yRho) {
        double h = Dataset_EVOL.getDepth(xRho, yRho, 0);
        return (h);
    }

// -----------------------------------------------------------------------------
    void auBord() {

        if (compteur_standstill > 100) {
            // System.out.println("MORT en STAND STILL");
            living = false;
        }
        if ((x > Dataset_EVOL.nx - 3.0f) || (x < 3.0f)) {
            living = false;
//             System.out.println("MORT AU BORD");
        }

        if ((y > Dataset_EVOL.ny - 3.0f) || (y < 3.0f)) {
            living = false;
//             System.out.println("MORT AU BORD");
        }

// TEST SUR LA PROFONDEUR :

        // 1 - Poisson dans la colonne d'eau (entre fond et surface)??
        if ((depth < bathyActuelle) || (depth > 0)) {
            living = false;
            // System.out.println("bathyActuelle = " + bathyActuelle + " ; depth = " + depth + "  -> MORT AU BORD");
            System.out.println("problem avec la profondeur de ce poisson. depth = " + depth);
            //System.out.println("bathyActuelle 2 = " + bathyActuelle);
        }
// 2 - Pronfondeur limite de la grille hydro?

        // AUTRES TEST :

        if (depth > 0) {
            living = false;
            System.out.println("depth = " + depth + "MORT EN L'AIR");
        }

        if (!Dataset_EVOL.isInWater(x, y)) {
//            System.out.println("MORT A TERRE");
            living = false;
        }

        if (Dataset_EVOL.isOnEdge(x, y)) {
            living = false;
            // System.out.println("Mort ON EDGE");
        }

        if (Dataset_EVOL.isCloseToCost(x, y)) {
            living = false;
            System.out.println("Mort CLOSE TO COAST");
        }
        // TEST DEBUG :
        // if (living == false) {
        //  System.out.println("Mort sur les bord de la grille : x = " + x + " , y = " + y + "bathyActuelle = " + bathyActuelle);
        //}
    }

    //--------------------------------------------------------------------------
    void stepSerial(int jour) {

// Indiv dans la zone consideree?
        //sinon on les y remet!
        if (y > Population.y_max) {
//              System.out.println("poisson quitte le domaine par le le Nord ");
//               y = Population.y_max - 1;
            Population.sortie_NORD++;

        }
        if (x > Population.x_max) {
//            x = Population.x_max-1;
//               System.out.println("poisson quitte le domaine par le l'Est ");
            Population.sortie_EST++;
        }
        if (y < Population.y_min) {
//            y = Population.y_min+1;
//               System.out.println("poisson quitte le domaine par le Sud ");
            Population.sortie_SUD++;
        }
        if (x < Population.x_min) {
//            x = Population.x_min+1;
//               System.out.println("poisson quitte le domaine par l'Ouest ");
            Population.sortie_OUEST++;
        }

// INCREMENT AGE -------------------------------------------------------
        age++;
///////////// ----- A U  B O R D  ? --------------------------------------------
// Elimine les individus se trouvant sur les bords du domaine :
        corr_depth();
        auBord();


        if (living) {

// DEPLACEMENT (Advection et Nage) ---------------------------------------------
            deplace_indiv(jour);

// MORTALITE PAR TEMPERATURE LETALE POUR LES OUEUFS ET LARVES :
            mortalite_ICHTYOPLANCTON_temperature_letale();

// PASSAGE EN COORD LON, LAT, DEPTH
            auBord();
            positGeog3D();
        }
    }

//------------------------------------------------------------------------------
    // POSITIONS POTENTIELLES DE PONTE :
    public void init_strates() {
        if (Simulation.Strategie.equals("HE")) {
            int d, i, j, k;

            strates =
                    new boolean[Simulation.nb_strat_temporelles][Population.X][Population.Y];

//    Nb_strates = 0;
            for (d = 0; d
                    < Simulation.nb_strat_temporelles; d++) {
                for (i = 0; i
                        < Population.X; i++) {
                    for (j = 0; j
                            < Population.Y; j++) {
//         System.out.println("i = " + i +" ; j = " + j + " ; d = " + d + " ; k = " + k);
                        strates[d][i][j] = false;
                    }
                }
            }
        }
    }

    public void erase_strates() {
        strates = null;

// Autre solution :
    /*
        public void add(Position_ponte p) throws IOException {
        this.positions_OK.add(p);
        }

        public void remove_positions() {
        this.positions_OK.clear();
        }
         */
    }

    void migration_ponte(int jour, int dt) {
        //  System.out.println("MIGRATION VERS ENVIRONNEMENT DE PONTE ");

        // 1 - Determiner l'env_ok le plus proche (zone de ponte potentielle la + proche)
        //      a) le même jour, dans un rayon de 10 km
        // en fait on cherche par resolution_temp donc le rayon = 10*resolution_temp
        int deplacement_max_par_jour_km = 10;
        int strate_temp = (int) Math.floor((jour) / ((float) Simulation.resolution_temporelle_scanEnv));

        float deplacement_max_dt = (float) deplacement_max_par_jour_km * ((float) dt / 24);
        // au moment "strate_temp" on cherche i et j tq
        //     strates[strate_temp][i][j]==true
        // &&  (i-Po.x)^2 + (j-Po.j)^2 < rayon_accessible_km*resolution_spatiale

        int rayon_recherche_ij = (int) Math.round(((double) deplacement_max_dt * ((float) Simulation.resolution_temporelle_scanEnv)) / (double) Population.mean_grid_size);
        int ii, i, jj, j;
        int i_target = 999;
        int j_target = 999;
        float dist_ij = 999;
        while (dist_ij == 999) {
            int ii_min = Math.max((int) x - rayon_recherche_ij, Population.i_min);
            int ii_max = Math.min((int) x + rayon_recherche_ij, Population.i_max);
            int jj_min = Math.max((int) y - rayon_recherche_ij, Population.j_min);
            int jj_max = Math.min((int) y + rayon_recherche_ij, Population.j_max);

            for (ii = ii_min + 1; ii < ii_max; ii++) {
                // PASSAGE AU COORDONEE DANS LE DOMAINE CHOISIT (de dimension X) :
                i = ii - Population.i_min;
                for (jj = jj_min + 1; jj < jj_max; jj++) {
                    // PASSAGE AU COORDONEE DANS LE DOMAINE CHOISIT (de dimension Y) :
                    j = jj - Population.j_min;

                    // System.out.println("jour = " + jour + " ; strate_temp" + strate_temp + " ; i = " + i + " ; j = " + j);
                    if (strates[strate_temp][i][j] == true) {
                        if (dist_ij > (x - i) * (x - i) + (y - j) * (y - j)) {
                            dist_ij = (float) Math.round((x - i) * (x - i) + (y - j) * (y - j));
                            i_target = ii;
                            j_target = jj;
                            System.out.println("Localisation d'un ENVIRONNEMENT DE PONTE a " + ((int) (double) dist_ij / (double) Population.mean_grid_size) + " km");
                        }
                    }
                }
            }
            // si on a pas trouver d'env_ok, alors on cherche plus loin, dans le prochain strate temporelle
            strate_temp++;
            if (strate_temp >= Simulation.nb_strat_temporelles) {
                living = false;
                break;
                // pas d'env de ponte, individu non selectionne.
                // (Il faudra considerer cette condition a partir d'un age = > juveniles,
                // car les juvenile ne font a priori pas de migration pour leur 1er ponte.)
            }
            rayon_recherche_ij += (int) Math.round(((double) deplacement_max_par_jour_km * ((float) Simulation.resolution_temporelle_scanEnv)) / (double) Population.mean_grid_size);
        }
        // 2 - nager dans cette direction SAUF si obstacle (côte)


        // ACODER -!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <---------------!!!
    }

    int nbre_au_hazard_entre_0_et_x(int x) {
        Random rand = new Random(); // constructeur
        int i = rand.nextInt(x + 1); // génération
        return (i);
    }

    void stepSerial_TEST_DEB(int jour) {

        //System.out.println(" TEST stepSerial_TEST_DEB S = " + S);
        CarCapa = 10000000; // Carying capacity en grammes : je met 10 tonnes
        int dt_jour = 1;


/// ----------------- P R E D A T I O N (selon la taille)-----------
        // OU COMPETITION
        // Inclue la mortalite par senescence
        // EN DEVELOPPEMENT 30 = taille max des larves..(?) NON c plutôt15-20 avec le modele de Vamara Koné
        // mortalité par "step_serial" (=jour)

        // voir script Matlab "Mortality.m" qui explore cette fonction
        // On re-initialise CarCapa :
        // Nouvelle Capacite de charge en bouffe du millieu
        CarCapa = 100000000; // grammes de Carbonne disponibles par km2 dans la couche des 0-20m

        //densite_biomasse_locale_par_km2 = la somme des biomasses des banc dans un rayon donne = maille de grille ROMS
        densite_biomasse_locale_par_km2_autres_super_individu = 0.f;
        densite_biomasse_de_ce_banc_par_km2 = (float) ((S * weight) / (Population.mean_grid_size * Population.mean_grid_size));
        densite_biomasse_locale_par_km2_TOTALE = densite_biomasse_locale_par_km2_autres_super_individu + densite_biomasse_de_ce_banc_par_km2;

        bathyActuelle = (float) getDepthBottom(x, y);

        float fish_length = DEB.getLength();
        if (fish_length < 1) {
            fish_length = 1;
        }
        S = Mortality.calcul_mortality(fish_length, DEB.getWeigth(), DEB.stade, age, CarCapa, S, densite_biomasse_locale_par_km2_TOTALE, densite_biomasse_de_ce_banc_par_km2, bathyActuelle, dt_jour);


        temp = 17;
        bouffe = 3 * DEB.X_K; // --> 100000 * DEB.X_K = bouffe non limitante (X_K = demisaturation)

        //            System.out.println(" TEST stepSerial_TEST_DEB temp =  " + temp + " bouffe = " +  bouffe + "  ((double) Simulation.dt_advec / 24) = " +  ((double) Simulation.dt_advec / 24));
        if (living) {

            // // ** OEUF ET LARVE :
            // : -------- A D V E C T I O N -------------------------------------
            if (age < Simulation.ageMinAtRecruitment) {

                //System.out.println("NEW DEPLACE ... jj = " + jj + "****************************");
                int it_debut, it_fin;
                it_debut =
                        (int) Population.nb_jours_ecoule_entre_roms_records * (24 * 60) / (int) Simulation.dt_advec;
                it_fin =
                        (int) (Population.nb_jours_ecoule_entre_roms_records + 1) * (24 * 60) / (int) Simulation.dt_advec;
                // *****************************************************
                int i = 0;
                for (int it = it_debut; it < it_fin; it++) {
                    double frac = it * Simulation.dt;
                    ///////////// ----- G R O W T H ------------------------------------------------
                    DEB.execute(temp, bouffe, (double) Simulation.dt_advec / (24 * 60));
                    Body_length = DEB.getLength(); // cm

//                    System.out.println(" Age = " + age + " jours et " + Simulation.dt_advec * i + " heures -> longueur = " + DEB.getLength() + " cm ");

                    i++;
                }
            }

            // ** JUVENILES ET ADULTE :
            // : -------- M O U V E M E N T S -------------------------------------
            // (UNIQUEMENT pour les individu qui ne sont pas en train de se reproduire)
            // --> ajouter un check d'etat de reproduction
            if (age >= Simulation.ageMinAtRecruitment && living) {
//System.out.println(" TEST stepSerial_TEST_DEB ADULTES");
                DEB.execute(temp, bouffe, 1); //<-- 1 car c'est le pas de temps de 1 jours
                Body_length = DEB.getLength(); // cm
            }
//            System.out.println(" Age = " + age + " jours -> longueur = " + DEB.getLength() + " cm  // Effectif = " + S + "  --> poids du Super individu = " + (S * DEB.getWeigth() / 1000000) + " tonnes ");

            // ** ADULTES MATURES :
            // : -------- R E P R O D U C T I O N -------------------------------------
            // 1 - Check si indiv mature et prêt à pondre
            // 2 - Si oui, PONDRE un nombre d'individu donné par DEB
            // + tard : check de l'environnement
            //          si envir_ok = true --> ponte
            //          si envir_ok = false --> MIGRATION VERS envir_ok
        }
        if (!living) {
            System.out.println(" STADE =  " + DEB.stade + "  -- MORT ");
        }
// INCREMENT AGE -------------------------------------------------------
        age++;
    }

    float CarCapa_locale() {

        // prend en entree la concentration de bouffe intstantanee et un pas de temps
        // donne en output la quantite de bouffe accessible a la predation pendant un temps dt
        // ----> IL FAUT PAR KM2 <---- pour pouvoir être comparé a densite_locale_par_km2

        // (en considerant que toute la bouffe est accessible dans un rayon de dt*longueur de poisson par seconde
        // ATTENTION cette approximation est de plus en plus fausse quand dt augmente (surestime le volume accessible)
        // dt EN MINUTES
//        float volume_explor_dt = (float) (Math.pow(Body_length, 3) * dt * 60) / 1000000; // Volume en m3 filtré par dt (en phase d'alimentation) = body_length^3 par seconde (approximation de la limite superieure;
        //       float CarCapa_locale = (densite_biomasse_plankton * volume_explor_dt); // en grammes
        //   float densite_biomasse_plankton = (float) (12 * densite_bouffe * 1000); // grammes par m3
        float CarCapa_locale;


        int carcapa_max_km2 = 100000000; // g soit 100 tonnes au km2

        CarCapa_locale = (float) (carcapa_max_km2 * f_reponse_fonctionnel);
        return CarCapa_locale;
    }

    void deplace_indiv(int jour) {
//    System.out.println("Deplacement des individus ");
        int threshold_metamorphose = 5;//
        float heure;
        double dt_frac; // 1/nbiter
        int dt_adapt = Simulation.dt_advec;
        double dx_swim = 0;
        double dy_swim = 0;
        double dx_tot, dy_tot;

        double dx_advec = 0;
        double dy_advec = 0;
//        double dz_advec = 0;

        float Body_length_TEST_KINESIS = 20;
        if (Simulation.TEST_KINESIS) {
            Body_length = Body_length_TEST_KINESIS;
        }

        boolean DT_ADAPTATIF = false;
        if (DT_ADAPTATIF) {
            if (Body_length > threshold_metamorphose) {
                dist_max_par_dt = 5;//Population.mean_grid_size; //km
                float vit_max_poisson = (float) ((Kinesis.Vitesse_max_Bls * Body_length / 100) * 3.6); // km/h
                dt_adapt = (int) (Math.round((dist_max_par_dt / vit_max_poisson) * 60f)); // minutes
            } else {
                dt_adapt = Simulation.dt_advec;
            }
        }
        float dt_jour = (float) ((float) dt_adapt / 60) / 24;

//// ------------ NE PAS TOUCHER C COMPLIQUE ET C VALIDE -----------------------
// INTERPOLATION DES DEPLACEMENT ENTRE LES SORTIES ROMS
// nb_jours_ecoule_entre_roms_records : dépend d'un compteur dans population...
        int it_debut, it_fin;
        it_debut = Population.nb_jours_ecoule_entre_roms_records * 24 * 60 / dt_adapt;
        it_fin = (Population.nb_jours_ecoule_entre_roms_records + 1) * 24 * 60 / dt_adapt;
//    System.out.println("it_debut = " + it_debut + " ; it_fin = " + it_fin);
        // *****************************************************
        dt_frac = 1.0f / Simulation.nbiter;

        for (int it = it_debut; it < it_fin; it++) {
            double frac = it * dt_frac;
// --FIN DU COMPLIQUE ET VALIDE (si,si) ----
            // Heure de la journee :
            float fraction_heure = ((float) dt_adapt) / 60;
            heure = it * fraction_heure + 1;
            while (heure > 24) {
                heure = heure - 24;
            }

// VERIFIONS QU'ON est dans l'eau : 
            if (!Dataset_EVOL.isInWater(x, y)) {
                System.out.println("MORT A TERRE 1111 BUGBUGBUGBUGBUG ");
            }
            try {
                // 1 - LECTURE/INTERPOLATION DES CHAMPS TEMPERATURE ET SALINITE :
                TS = Dataset_EVOL.getFields_SaltTemp(x, y, z, frac);
                grad_temp = (TS[0] - temp) / dist_parcourue_km;
                //System.out.println("dist_parcourue_km = " + dist_parcourue_km + " ; grad_temp = " + grad_temp);

                if (isNaN(grad_temp)) {
                    grad_temp = 0.001;
                }
                temp = TS[0];
                salt = TS[1];

                if (age == 1) {
                    imprinting();
                }


                // 2 - LECTURE/INTERPOLATION DES CHAMPS PLANKTON ET OXYGENE :
//                int nz = Dataset_EVOL.get_nz(); //(nz : nbre de niveau sigma, pour avoir la temperature de surface)
                get_bouffe(frac); // actualise "bouffe"

                // LECTURE/INTERPOLATION DES COURANTS
                // SI ADVECTION = ON : LECTURE/INTERPOLATION DES CHAMPS DE VITESSE U, V et W :
                if (Simulation.ADVECTION) {
                    uv = Dataset_EVOL.getFields_uv(x, y, z, frac, dt_adapt);
                    dx_advec = uv[0];
                    dy_advec = uv[1];

                    if ((dx_advec == 0) && (dy_advec == 0)) {
                        //System.out.println("dx_advec = " + dx_advec + " ; dy_advec = " + dy_advec);
                        compteur_standstill = compteur_standstill + 1;
                    }
//                    dz_advec = uvw[2];

                    /* TEST                   if(frac == it_debut*dt_frac){
                    System.out.println("Jour = " + jour + " - dx_advec = " + dx_advec);
                    System.out.println("Jour = " + jour + " - dy_advec = " + dy_advec);
                    }*/
                }
            } catch (IOException e) {
                System.out.println("youstone on a un probleme dans deplace_individu : " + e.getMessage());
            }


/// ----------------- P R E D A T I O N (selon la taille)-----------
            // OU COMPETITION
            // Inclue la mortalite par senescence
            // EN DEVELOPPEMENT 30 = taille max des larves..(?) NON c plutôt15-20 avec le modele de Vamara Koné
            // mortalité par "step_serial" (=jour)

            // voir script Matlab "Mortality.m" qui explore cette fonction
            // On re-initialise CarCapa :
            // Nouvelle Capacite de charge en bouffe du millieu
            CarCapa = CarCapa_locale(); // grammes de Carbonne disponibles par km2 dans la couche des 0-20m

            //densite_biomasse_locale_par_km2 = la somme des biomasses des banc dans un rayon donne = maille de grille ROMS
            densite_biomasse_locale_par_km2_autres_super_individu = (float) ((Population.Biom_tot_2D[(int) x][(int) y]) / (Population.mean_grid_size * Population.mean_grid_size));
            densite_biomasse_de_ce_banc_par_km2 = (float) ((S * weight) / (Population.mean_grid_size * Population.mean_grid_size));
            densite_biomasse_locale_par_km2_TOTALE = densite_biomasse_locale_par_km2_autres_super_individu + densite_biomasse_de_ce_banc_par_km2;

            bathyActuelle = (float) getDepthBottom(x, y);

            float fish_length = DEB.getLength();
            if (fish_length < 1) {
                fish_length = 1;
            }
            S = Mortality.calcul_mortality(fish_length, DEB.getWeigth(), DEB.stade, age, CarCapa, S, densite_biomasse_locale_par_km2_TOTALE, densite_biomasse_de_ce_banc_par_km2, bathyActuelle, dt_jour);
            Population.Biom_landings_2D[(int) x][(int) y] = Mortality.Biom_landings;
            
            // SUPRESSION DU SUPER INDIVIDU SI S < 10 ---------------------------------------
            if (S < 1000) { // car en dessous de 1000 on estimme que ça devient compliqué au niveau de la prédation (pour former des bancs)
                living = false;
                //        System.out.println("DISPARITION DU BANC DE POISSON num" + id + " --> CarCapa = " + CarCapa + " ;  mortalite predation par jour = " + Mortality.M_predation + " ;  mortalite senescence par jour = " + Mortality.M_senescence);
            }

            //         System.out.println("Body_length = "+ Body_length+ " , S = " + S + " MORTALITE = " +  mort + " --> new S = " + (int) Math.round(S * Math.exp(-mort)) );
            //        S = (int) Math.round(S * Math.exp(-mort)); (WATKINS ET ROSE)
            //        S = (int) Math.round(S * (1-mort)); //(Pierre Auger)

// ---------------------- PECHE (selon taille et selon zone) -------------------
            // A CODER --
// -------------------FIN PREDATION --------------------------------------------

// CALCULE DE LA NAGE HORIZONTALE PAR EXTENTYED KINESIS
            if (Simulation.TEST_KINESIS == true) {
                //depth = -1; // on tester le suivi du champs de bouffe de surface
                temp = 20; // Pour que le champs de SST n'interfere pas : la qualité de l'habitat ne dépend que de la bouffe
                Body_length = Body_length_TEST_KINESIS; // consideron la nage de poissons de 20cm.
                do_kinesis(dt_adapt);
                Vertical_behavior(heure);
            } else if (Body_length > threshold_metamorphose) {
                if (DEB.adult_ready_to_spawn && Simulation.Strategie.equals("HE") && (first_spawn == false)) {
                    migration_ponte(jour, dt_adapt);
                    // A CODER POUR FOURNIR dx_kinesis_swim_km et dy_kinesis_swim_km
                } else {
                    do_kinesis(dt_adapt); // afecte des valeure a dx_kinesis_swim_km et dy_kinesis_swim_km
                }
            }
            // ) distance parcourue durant le pas de temps, en nombre de maille de la grille :
            dx_swim = (double) dx_kinesis_swim_km / (double) Population.mean_grid_size;
            dy_swim = (double) dy_kinesis_swim_km / (double) Population.mean_grid_size;

// CALCULE DE LA NAGE HORIZONTALE PAR LUTTE CONTRE LE COURANT
// C'est un pourcentage de la nage totale qui doit s'orienter contre le courant


            //--------- DEPLACEMENT HORIZONTAL DU POISSON : NAGE + AVECTION ----------------
// ajoutons un facteur de lutte contre le courant de 1 Bls
//

            /*            float nage_cc = 0;
            if (Body_length > threshold_metamorphose){
            nage_cc = 1.f; // 0 = Resistance au courant
            // A FAIRE VARIER + SELON LA TAILLE DU POISSON??
            }
            float dx_advec_cor = (float) (nage_cc*dx_advec);
            float dy_advec_cor = (float) (nage_cc*dy_advec);           
            dx_tot = dx_advec_cor + dx_swim; // ; //+
            dy_tot =  dy_advec_cor + dy_swim; //dy_swim +
             */
            dx_tot = dx_advec + dx_swim; // ; //+
            dy_tot = dy_advec + dy_swim; //dy_swim +
// Verification que le mouvement de NAGE + ADVECTION ne nous amene pas à terre :
            if (Dataset_EVOL.isInWater(x + dx_tot, y + dy_tot) //) {
                    && (!Dataset_EVOL.isCloseToCost(x + dx_tot, y + dy_tot))) {
//                   && (!Dataset_EVOL.isOnEdge(x + dx_tot, y + dy_tot))) {
                // dans ce cas on effectue le deplacement
                x = x + dx_tot;
                y = y + dy_tot;
                dist_parcourue_km = (float) Math.sqrt(dx_tot * dx_tot + dy_tot * dy_tot) * Population.mean_grid_size;
                if (dist_parcourue_km > 0) {
                    compteur_standstill = 0;
                }

            } else {
                // sinon on s'arrette (stand still) et on met a zero la vitesse de kinesis
                V_kinesis[0] = 0;
                V_kinesis[1] = 0;
                dist_parcourue_km = 0;
                compteur_standstill = compteur_standstill + 1;
                //
            }
// On reverifie qu'on est toujours dans l'eau =
            if (!Dataset_EVOL.isInWater(x, y)) {
                System.out.println("MORT A TERRE 2222 BUGBUGBUGBUGBUG ");
            }
// RETEST QUE LA KINESIS NOUS A PAS AMENER AU BORD ----------------
            auBord();
            if (living == false) {
                break;
            }

//--------- DEPLACEMENT VERTICAL DU POISSON : NAGE + AVECTION ----------------

            bathyActuelle = (float) getDepthBottom(x, y);

//            if (z + dz_advec < Dataset_EVOL.nz) {
            // dans ce cas on effectue le deplacement
            //              z = z + dz_advec;
            //        }
            // On néglige l'advection VERTICALE

            if (age <= Simulation.egg_duration) {
                Egg_buoyancy();
            } else if (Body_length < threshold_metamorphose) {
                DVM_larves();
            } else {
                Vertical_behavior(heure);
            }
            // Correctif sur la profondeur
            corr_depth();

// CROISSANCE ----- G R O W T H ------------------------------------------------
// (SI PB de calcul moyenner temp et bouffe sur 24 heures et appeler DEB.execute 1 fois par jour dans step.serial)
            double dt_en_jours = ((double) dt_adapt / 24) / 60;
            DEB.execute(temp, bouffe, dt_en_jours);
            Body_length = DEB.getLength(); // cm
        }
    }

    void Egg_buoyancy() {
///////////// ----- E G G  B U O Y A N C Y -------------------------------------
        if (Simulation.flagBuoy) {
            //System.out.println("util.BuoyancyScheme.waterDensity(salt_init, temp_init) : " + util.BuoyancyScheme.waterDensity(salt_init, temp_init));
            // util.BuoyancyScheme.move : Nous donne le mouvement en m pour un pas de temps de Simulation.dt_advec
            double buoy_move = util.BuoyancyScheme.move(egg_density, salt, temp);
            //System.out.println("buoy_move en metres : " + buoy_move);
            depth = (float) (Math.min(-0.1f, depth + buoy_move)); // buoyancyMeters is in m/ss

            if (depth < bathyActuelle) {
                depth = bathyActuelle + 1;
            }
            //System.out.println("depth 2 : " + depth);
        }
    }

    void DVM_larves() {

///////////// ----- S W I M ( D V M ) ------------------------------------------
// MIGRATION VERTICALE (DVM SURFACE - Oxycline) :
        if (Simulation.flagSwim_random) {
            depth = (float) (-Math.random() * Simulation.Swim_lowerlimit);
        }
        if (depth < bathyActuelle) {
            depth = bathyActuelle + 1;
        }
    }

///////////// ----- S W I M  KINESIS ***   -------------------------------------
    void do_kinesis(int dt_adapt) {
        //System.out.println("On nage kinesis : Body_length = " + Body_length + " , CarCapa = " + CarCapa + " , densite_biomasse_locale_par_km2 = " + densite_biomasse_locale_par_km2);

        // (bouffe_max_0 + bouffe_max_1 + bouffe_max_2 + bouffe_max_3));
        // bouffe max : pour faire comme si yavait une memoire des endroit les mlieux qui puissent exister et que ça change le comportement..
        // APRES ON LA CALCULERA SELON LES INDICATIONS DE PEDRO

        //weight_old = weight;
        //weight = DEB.getWeigth();
        //stade = DEB.stade;

        /*       if (weight-weight_old <0){
        System.out.println("ATTENTION : POISSON QUI MAIGRIT! ");
        }
         */
        double T_opt = temperature_natal;

        // TOLERANCE AUGMENTE AVEC AGE :     
        //double sigma_t = Math.log(Body_length / 3);

        // TOLERANCE DIMINUE AVEC AGE : 
        double sigma_t = Math.min(4 - Math.log(Body_length / 3), 0.5);
       

        //        double sigma_t = Math.log(age/30+1);
        // age de 1 mois --> sigma_t = 0.7
        // age de 3 mois --> sigma_t = 1.4
        // age de 6 mois --> sigma_t = 1.9
        // age de 1 ans --> sigma_t = 2.5
        // age de 2 ans --> sigma_t = 3.2
        // age de 3 ans --> sigma_t = 3.6
        // age de 4 ans --> sigma_t = 3.9
        // age de 5 ans --> sigma_t = 4.1

        Qold = Q;
        Q = Kinesis.Habitat_quality(densite_biomasse_locale_par_km2_autres_super_individu,
                CarCapa, Body_length, temp, grad_temp, f_reponse_fonctionnel, T_opt, sigma_t,
                lon, lat, (int) bathyActuelle);

        V_kinesis_old = V_kinesis;
        V_kinesis = Kinesis.move(V_kinesis_old, Body_length, Q, Qold);

//          System.out.println("V_kinesis (x,y,z) : " + V_kinesis[0] + " , " + V_kinesis[1] + " , " +V_kinesis[2]);
// ON ACTUALISE LA POSITION :
// pas de temps = dt_advec en heures (pour l'instant)
// V_kinesis en m/s
// il faut convertir le deplacement Kinesis de metres en coord grille 

        // a) distance parcourue durant le pas de temps, en km :
        float dt_en_heure = ((float) dt_adapt) / 60;

        dx_kinesis_swim_km = V_kinesis[0] * 3.6 * dt_en_heure;
        dy_kinesis_swim_km = V_kinesis[1] * 3.6 * dt_en_heure;
        // dz_swim_m = V_kinesis[2] * 3600 * dt_en_heure/100;


        //       System.out.println("Poisson de " + Body_length + "cm ,  Q =  " + Q + " , vitesse_x = "  + Math.round(V_kinesis[0]*3.6) +"km/h");
//            System.out.println("Poisson de " + Body_length + "cm - vitesse_x = "  +V_kinesis[0]*3.6 +"km/h  -  dx_kinesis_swim_km =  " + dx_kinesis_swim_km + " km en " + Simulation.dt_advec + "heure soit " + Vx_kinesis_grd + " point de grille");
//            System.out.println("vitesse_y = "  +V_kinesis[1]*3.6 +"km/h  -  dy_kinesis_swim_km =  " + dy_kinesis_swim_km + " km en " + Simulation.dt_advec + "heure soit " + Vy_kinesis_grd + " point de grille");
//            System.out.println("dz_swim_m =  " + dz_swim_m + " m en " + Simulation.dt_advec + "heure");

    }

    void Vertical_behavior(float heure) {
// COMPORTEMENT VERTICAL DES BANCS : depend de l'heure et de la bathy
// PROF MAXIMALE OU SONT OBSERVES DES BANCS? Smith et Brown 2002 : quasiment rien en dessous de 1000m..
        // ATTENTION : dephth et bathyActuelle = valeures NEGATIVE!!!

        bathyActuelle = (float) getDepthBottom(x, y); //Dataset_EVOL.getDepth(x, y, 0);

        float pic_prof;
        float variance_prof = 10;

        if (heure <= 12) { // 12h de jour : bancs + en profondeur
            pic_prof = 30;//12;
        } else { // 12h de nuit : bancs + en surface
            pic_prof = 20;//6;
        }

        double prof_abs = Math.abs(Kinesis.getGaussian(pic_prof, variance_prof));
        while ((Math.abs(bathyActuelle) - prof_abs < 0)) {
            prof_abs = Math.abs(2 * Math.abs(bathyActuelle) - prof_abs);
        }
        depth = (float) -prof_abs;
//for(int i=1;i<100;i++){
//double prof_test = Math.abs(Kinesis.getGaussian(pic_prof, variance_prof));
//System.out.println( " prof_abs =  " + prof_abs + " , bathyActuelle = " + bathyActuelle);
//}
        //  System.out.println(" Heure =  " + heure + " ; bathy = " +bathyActuelle + " ; prof du banc = " + depth);
// TEST
        if ((depth < bathyActuelle) || (depth > 0)) {
            living = false;
            // System.out.println("bathyActuelle = " + bathyActuelle + " ; depth = " + depth + "  -> MORT AU BORD");
            System.out.println(" TEST 1 problem avec la profondeur de ce poisson. depth = " + depth);
            //System.out.println("bathyActuelle 2 = " + bathyActuelle);
        }
//FIN TEST

        /*
        try {
        positGeog3D(); // Passage des coordonées grille aux coordonnée lon lat depth
        } catch (Exception e) {
        System.out.println("probleme x= " + x + " , y= " + y + "z= " + z);
        System.out.println("probleme positGeog3D : " + e.getMessage());
        }
         *
         */
    }

    void corr_depth() {
// DEBUT Correction sur les valeurs aberrantes de depth : ----------------------
        bathyActuelle = (float) getDepthBottom(x, y); //Dataset_EVOL.getDepth(x, y, 0);

        //   ss = 0;
        //   while ((depth < bathyActuelle + 1) || (depth > -0.1)) {
        //      ss++;

        if (depth > -0.1) {
            depth = -1.f;//(float) (-0.1 - Math.random() * 3);
        } else if (depth < bathyActuelle + 1) {
            depth = bathyActuelle + 2.f; //(float) (bathyActuelle + 2 + Math.random() * 3);
        }

        //    if (ss > 100) {
        //      System.out.println(" impossible (corr_depth) - bathyActuelle = " + bathyActuelle);
        //    System.out.println(" x = " + x + " y = " + y + " depth = " + depth);
        //  System.exit(0);
        //}
        //  }
/// FIN de correction sur les valeurs aberrantes de depth ----------------------}
// Passage de la profondeur en coordonnées Grille ******************************
        z = Dataset_EVOL.depth2z(x, y, depth);
    }

    void get_bouffe(double frac) {
        double[] pGrid = new double[3];
        pGrid[0] = x;
        pGrid[1] = y;
        pGrid[2] = 31;//z;


        if (Simulation.TEST_KINESIS) {
            pGrid[2] = 31;
            //          pGrid[2] = Dataset_EVOL.depth2z(x, y, 0);;//31; // on vuet tester le suivi du champs de bouffe de surface
//double [] test_plankton_surf = Dataset_EVOL.getPlankton_PISCES(pGrid, frac);
//System.out.println("test_plankton_surf " + test_plankton_surf[1] + " z = " + pGrid[2]);
            //        pGrid[2] = Dataset_EVOL.depth2z(x, y, bathyActuelle );;//31; // on vuet tester le suivi du champs de bouffe de surface
            //          double [] test_plankton_fond = Dataset_EVOL.getPlankton_PISCES(pGrid, frac);
//System.out.println("test_plankton_fond " + test_plankton_fond[1] + " z = " + pGrid[2]);
        }

        plankton = Dataset_EVOL.getPlankton_PISCES(pGrid, frac);
        // PISCES : {Diatoms, NanoPhyto, MicroZoo, MesoZoo, O2}
        bouffe = plankton[0] + plankton[1] + plankton[2] + plankton[3];

        f_reponse_fonctionnel = Math.max(0, bouffe / (bouffe + DEB.X_K));
//        if (f_reponse_fonctionnel<0){
//          System.out.println("BUG f_reponse_fonctionnel = " + f_reponse_fonctionnel );
//        }

        /*
        if (Simulation.TEST_KINESIS){
        if (lat>20){
        //System.out.println("ZONE RICHE ARTIFICIELLE " );
        bouffe = 10*DEB.X_K; // pour test accumulation des poissons au sud de 18N?
        }
        else {bouffe = 0.5*DEB.X_K;}

        }
        /* // POUR si on veut garder la memoire des concentration de bouffe les plus importantes rencontrées :
        if (plankton[0] > bouffe_max_0) {
        bouffe_max_0 = plankton[0];
        }
        if (plankton[1] > bouffe_max_1) {
        bouffe_max_1 = plankton[1];
        }
        if (plankton[2] > bouffe_max_2) {
        bouffe_max_2 = plankton[2];
        }
        if (plankton[3] > bouffe_max_3) {
        bouffe_max_3 = plankton[3];
        }
         */
    }

    void imprinting() {

        //*******************************************
        // environnement initial : (POUR HOMING ENVIRONEMENTAL)
        temperature_natal =
                TS[0];
        salinite_natal =
                TS[1];
        profThermo_natal =
                TS[2];
        phyto_init =
                0; // reste a rajouter
        egg_density = util.BuoyancyScheme.waterDensity(salinite_natal, temperature_natal) - Simulation.egg_excess_buoyancy;
        //System.out.println("util.BuoyancyScheme.waterDensity(salt_init, temp_init) =  " + util.BuoyancyScheme.waterDensity(salt_init, temp_init));
        //System.out.println("egg_density = " + egg_density);

        //           System.out.println("-----> Imprinting finis ");
    }

    void mortalite_ICHTYOPLANCTON_temperature_letale() {

        if (age < 30) {
            ///////////// ----- L E T H A L  T E M P E R A T U R E -------------------------
            //                System.out.println("temp : " + temp);
            isdead_temp = ((temp < Simulation.temp_lethal_min || //PAR TEMPERATURE MIN
                    temp > Simulation.temp_lethal_max));    //PAR TEMPERATURE MAX

            // Juste pour enregistrer les temperatures extremes rencontrées :
            if (tempmin > temp) {
                tempmin = (float) temp;
            }
            if (tempmax < temp) {
                tempmax = (float) temp;
            }

            if (isdead_temp) {
                living = false;
            }
        }
    }

    void get_identity() {
/*
 * Idee de mettre des info sur le lieu de naissance dans l'id des poissons..:
        int lat_id = (int) Math.round(lat_init * 10);
//    int lon_id = (int) Math.abs(Math.round(lon_init*10));
        if (day_init < 10) {
            id = Long.valueOf("" + lat_id + "00" + day_init + Population.compteur_id_poissons);
        } else if (day_init < 100) {
            id = Long.valueOf("" + lat_id + "0" + day_init + Population.compteur_id_poissons);
        } else {
            id = Long.valueOf("" + lat_id + day_init + Population.compteur_id_poissons);
        }
 *
 */
        id =Population.compteur_id_poissons;
        Population.compteur_id_poissons++;
    }
}
