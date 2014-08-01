package evolution;

/**
 *
 * @author timbrochier
 */
public class Senegal extends Dataset_EVOL {


// [EVOL2010] : Les parametres sont presque tous lu dans le fichier config.
// Seul reste a specifier ici :
//        - la liste des annees et celle des profondeur potentielles (pour
    // homing environnemental) car pas trouver comment lire une liste dans le
// // fichier config
//        - le Nom des variables dans le netcdf output de roms, car ca ne change
    // presque jamais (getFieldsName)
    public Senegal() {

        // INITIALISATION DU FICHIER CONFIG :
        util.Lire_configuration param_config = new util.Lire_configuration();
        try {
            String fichier_config = System.getProperty("user.dir") + "/" + IBM.config_file;
            util.Lire_configuration.Setup(fichier_config);
        } catch (Exception e) {
// TODO Auto-generated catch block
            System.out.println("probleme dans initialisation du fichier config : " + e.getMessage());
            e.printStackTrace();
        }

        // Lecture du fichier config :
        try {
            //region = "PERU_actuel_moyenne";
            // PERU_PISCES_sixieme , PERU_sixieme , PERU_neuvieme, PERU_preindus, PERU_4xCO2, PERU_actuel, PERU (ex PERU_INTERVINCENT_CLIMPIERRICK)
            // Pour la grille roms Peru10,  1 maille =9km (1/9°)
            // Pool de sorties hydro : (A) 4 ans de Climato Pierrick, year 12-15
            //                         (B) 8 ans de Interanuel Vincent, 1992-2000

            region = param_config.getConfig_string("region");
            //System.out.println ("region = " + region);

            SIMU = param_config.getConfig_int("SIMU");

            ////CONDITIONS INITIALES //////
            // ZONE DE PONTE PERU:
            // Cette zone doit largement englober la zone de ponte
            // --> ou alors rajouter automatiquement un "sponge" dans Population
            lat_min = param_config.getConfig_double("lat_min");
            lat_max = param_config.getConfig_double("lat_max");
            lon_min = param_config.getConfig_double("lon_min");
            lon_max = param_config.getConfig_double("lon_max");
            bathy_max = param_config.getConfig_int("bathy_max");

            prof_ponte_min = param_config.getConfig_int("prof_ponte_min"); //  PROFONDEUR DE PONTE INITIALE MAX
            prof_ponte_max = param_config.getConfig_int("prof_ponte_max");//  PROFONDEUR DE PONTE INITIALE MIN

            // Profondeur du talu continental (pour le critere de retention)
            prof_talu = param_config.getConfig_int("prof_talu"); // [metres]

            // SPONGE ou TAMPON : distance en km entourant la zone de ponte, tel que les individu ne puissent pas en sortir
            sponge_km = param_config.getConfig_int("sponge_km"); // Km

            // PATH DES SORTIES HYDRO ROMS:
            directory_roms = param_config.getConfig_string("directory_roms") + region + "/";

            //directory_roms = "/Users/timbrochier/Documents/Sorties_Models/" + region + "/";
            // PECH12 :
            //directory_roms = "/Users/timbrochier/Documents/Sorties_Models/Pech12/";

            directory_Suplementary_data = param_config.getConfig_string("directory_Suplementary_data");
            main_output_directory = param_config.getConfig_string("main_output_directory");


            PonteOBS_filename = param_config.getConfig_string("PonteOBS_filename");
            PonteCLM_filename = param_config.getConfig_string("PonteCLM_filename");
            oxyclin_filename = param_config.getConfig_string("oxyclin_filename");


            Clim_chlaSeaWiFS = param_config.getConfig_string("Clim_chlaSeaWiFS");

            sufixe = param_config.getConfig_string("sufixe");

        } catch (Exception e) {
// TODO Auto-generated catch block
            System.out.println("probleme dans Lire_configuration : " + e.getMessage());
            e.printStackTrace();
        }


        // LISTE YEAR :

        if (region.equals("SENEGAL")) {
            // liste PERU_PISCES (jeudi 4 juin 2009)
//            yearlist100 = new int[]{6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
//            yearlist100 = new int[]{1980, 1981,1982,1980, 1981,1982,1980, 1981,1982,1980, 1981,1982,1980, 1981,1982,1980, 1981,1982};
            yearlist100 = new int[]{1980, 1981, 1982,1983,1984, 1985,1986,1987, 1988,1989,1990, 1991,1992,1993, 1994,1995,1996, 1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009};

        } else if (region.equals("PERU")) {
            yearlist100 = new int[]{1995, 12, 13, 15, 1999, 1992, 1994, 1993, 1992, 1997, 15, 1993, 15, 1992, 1998, 1994, 1996, 13, 12, 1996, 1996, 13, 1992, 15, 1994, 1994, 1999, 14, 13, 14, 1994, 1999, 14, 1996, 2000, 1998, 1999, 1995, 14, 1992, 2000, 12, 1993, 12, 1999, 2000, 1992, 1994, 1997, 1992, 1997, 14, 15, 1995, 1992, 1997, 12, 1998, 1998, 1997, 12, 1999, 15, 13, 1996, 12, 1996, 1995, 2000, 1993, 1997, 1993, 1998, 2000, 15, 1994, 14, 1997, 1995, 1993, 1993, 14, 1994, 12, 1993, 13, 13, 1999, 2000, 1998, 14, 1996, 13, 1995, 1995, 1996, 1995, 1998, 15, 2000};
        }

        // POUR OPPORTUNISME  et Homing Environnemental :
        //   ZONE DE PONTE POTENTIELLE: zone initiale
        //   PROFONDEUR DES ZONES DE PONTE POTENTIELLE :
//        prof_potentielles = new int[]{-5, -10, -15, -20, -25, -30, -35, -40, -45, -50};
        //prof_potentielles = new int[]{-5, -15};
        prof_potentielles = new int[]{-10};
    }

    void getFieldsName() {

        // Nom des variables dans le netcdf output de roms:
        strXiDim = "xi_rho";
        strEtaDim = "eta_rho";
        strZDim = "s_rho";
        strTimeDim = "time";
        strLon = "lon_rho";
        strLat = "lat_rho";
        strBathy = "h";
        strMask = "mask_rho";
        strU = "u";
        strV = "v";
        strW = "w";
        strOmega = "omega";
        strZeta = "zeta";
        strTp = "temp"; //scrum_time
        strSal = "salt";
        strTime = "time";
        strKv = "AKt";
        strHBL = "hbl";

        strPn = "pn";
        strPm = "pm";
        strThetaS = "theta_s";
        strThetaB = "theta_b";
        strHc = "hc";

// PISCES :

        strDiatoms = "DIA";
        strNanoPhyto = "NANO";
        strMicroZoo = "ZOO";
        strMesoZoo = "MESO";
        strO2 = "O2";
/*
// NPZD :
        strPhyto = "PHYTO";
        strZoo = "ZOO";
        strDet = "DET";
        strNO3 = "NO3";
        strCHLA = "CHLA";

        strDiatomsChl = "DCHL";
        strNanoPhytoChl = "NCHL";
*/
    }

}
