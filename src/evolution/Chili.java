package evolution;

public class Chili extends Dataset_EVOL {

    public Chili() {

        region = "CHILI";
        // Pour la grille roms Peru10,  1 maille =9km (1/9°)
        // Pool de sorties hydro : (A) 4 ans de Climato Pierrick, year 12-15
        //                         (B) 8 ans de Interanuel Vincent, 1992-2000

        // PATH DES SORTIES HYDRO ROMS:
        directory_roms = "/Volumes/SilverTouch/ROMS/PERU_INTERVINCENT_CLIMPIERRICK/";
        sufixe = ".nc";

        // METTRE ICI LE NUMERO DE LA SIMULATION
        SIMU = 1;

        ////CONDITIONS INITIALES //////
        // ZONE DE PONTE PERU: 
        lat_min = -22;
        lat_max = -4;
        lon_min = -85;
        lon_max = -70;
        bathy_max = 9999;

        prof_ponte_min = 0; //  PROFONDEUR DE PONTE INITIALE MAX 
        prof_ponte_max = 30;//  PROFONDEUR DE PONTE INITIALE MIN

        // Profondeur du talu continental (pour le critere de retention)
        prof_talu = 500; // [metres]

        // Liste de 100 ans : tirage aleatoire dans le pool des années hydro disponible (avec matlab) 
        //Liste avec Clim + 1992 a 2000 
        yearlist100 = new int[]{1995,12,13,15,1999,1992,1994,1993,1992,1997,15,1993,15,1992,1998,1994,1996,13,12,1996,1996,13,1992,15,1994,1994,1999,14,13,14,1994,1999,14,1996,2000,1998,1999,1995,14,1992,2000,12,1993,12,1999,2000,1992,1994,1997,1992,1997,14,15,1995,1992,1997,12,1998,1998,1997,12,1999,15,13,1996,12,1996,1995,2000,1993,1997,1993,1998,2000,15,1994,14,1997,1995,1993,1993,14,1994,12,1993,13,13,1999,2000,1998,14,1996,13,1995,1995,1996,1995,1998,15,2000};
        
        // FRANCOIS, CI DESSOUS NE T'EN OCCUPE PAS (PAS POUR LE NATAL HOMING)

        // POUR OPPORTUNISME  et Homing Environnemental :
        //   ZONE DE PONTE POTENTIELLE: zone initiale
        //   PROFONDEUR DES ZONES DE PONTE POTENTIELLE : 
        prof_potentielles = new int[]{-5, -10, -15, -20, -25, -30, -35, -40, -45, -50, -55};
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
        strOmega = "omega";
        strZeta = "zeta";
        strTp = "temp";
        strSal = "salt";
        strTime = "scrum_time";
        strKv = "AKt";
        strHBL = "hbl";

        strPn = "pn";
        strPm = "pm";
        strThetaS = "theta_s";
        strThetaB = "theta_b";
        strHc = "hc";

    //      strLargePhyto = GrowthModel.LARGE_PHYTO;
    //      strLargeZoo = GrowthModel.LARGE_ZOO;
    //      strSmallZoo = GrowthModel.SMALL_ZOO;

    }
    //---------- End of class
}

