package util;

/** import java.io */
//import java.io.File;
import java.io.IOException;
import ucar.ma2.Index;

/** import java.util */
//import java.util.ArrayList;

/** import netcdf */
import ucar.ma2.Array;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.dataset.NetcdfDataset;
import evolution.Dataset_EVOL;

/**
 *
 * @author timbrochier
 */
public class Get_SpawnFile {

    static NetcdfFile nc_OBS, nc_CLM;
    static int NB_YEAR,  nx_obs,  ny_obs;
    static int[] liste_year_obs;
    static int year_debut,  i_deb;
    static int year_fin,  i_fin;
    static double[] lon_obs;
    static double[] lat_obs;
    static int[] month_obs;
    static float[][][][] eggs_obs;
    static float[][][] eggs_CLM;

    static int ilat_deb,  ilat_fin;
    static double[] lon_ponte;
    static double[] lat_ponte;
    static double[] densite_ponte;
    static int n_ponte;
    public static double[] lon_ponte_3mois;
    public static double[] lat_ponte_3mois;
    public static double[] densite_ponte_3mois;

    public static int n_ponte_3mois;

    public static void init() throws IOException {

        //long t0 = Simulation.get_t0();
        // ici on ne repere plus les fichiers par leur "t0" mais par leur nom directement:
        // open(getFile(t0));
        String PATH_SpawnFile = Dataset_EVOL.directory_Suplementary_data + Dataset_EVOL.PonteOBS_filename;

        String PATH_SpawnFile_CLIMATO = Dataset_EVOL.directory_Suplementary_data + Dataset_EVOL.PonteCLM_filename;


        System.out.println("lecture des Obs + CLIMATO... ");

        System.out.println("PATH_SpawnFile : " + PATH_SpawnFile);
        System.out.println("PATH_SpawnFile_CLIMATO : " + PATH_SpawnFile_CLIMATO);

        nc_OBS = NetcdfDataset.openFile(PATH_SpawnFile, null);
        nc_CLM = NetcdfDataset.openFile(PATH_SpawnFile_CLIMATO, null);



//        nbTimeRecords = ncIn.findDimension(strTimeDim).getLength();

        lecture_Obs();
        lecture_CLM();
    }

    /**
     *INIFile =  file contenant les dates et lieu de d observation de la ponte.

     * Reads the properties under the RELEASE section, when the release mode
     * from NetCDF file is selected.
     *
     * @param file an INIFile, the configuration file.
     * @param section a String, the name of the section
     */
    static void lecture_Obs() throws IOException {

        int[] shape;
        int[] origin;

        int latlim_nord = -5;
        int latlim_sud = -15;

        Array arrYear, arrLon, arrLat;
        Index index;

        NB_YEAR = 0; // number of year with observations

        try {
            NB_YEAR = nc_OBS.findDimension("year").
                    getLength();
            nx_obs = nc_OBS.findDimension("longitude").
                    getLength();
            ny_obs = nc_OBS.findDimension("latitude").
                    getLength();

            System.out.println("NB_YEAR = " + NB_YEAR + " ; nx_obs = " + nx_obs + " ; ny_obs = " + ny_obs);


// GET LIST YEAR :
            ucar.nc2.Variable y = nc_OBS.findVariable("year");
            shape = y.getShape();      // get copy of shape array
            origin = new int[y.getRank()]; // start with all zeroes
            arrYear = y.read(origin, shape); //.reduce().copyToNDJavaArray();    // read data from file
            liste_year_obs = new int[NB_YEAR];
            index = arrYear.getIndex();
            for (int j = 0; j < NB_YEAR; j++) {
                index.set(j);
                liste_year_obs[j] = arrYear.getInt(index);
            //   System.out.println("liste_year_obs[j] = " + liste_year_obs[j]);
            }


// GET LIST LON : 
            ucar.nc2.Variable lo = nc_OBS.findVariable("longitude");
            shape = lo.getShape();      // get copy of shape array
            origin = new int[lo.getRank()]; // start with all zeroes
            arrLon = lo.read(origin, shape);
            lon_obs = new double[nx_obs];
            index = arrLon.getIndex();
            for (int j = 0; j < nx_obs; j++) {
                index.set(j);
                lon_obs[j] = arrLon.getDouble(index);
            }

// GET LIST LAT :
            ucar.nc2.Variable la = nc_OBS.findVariable("latitude");
            shape = la.getShape();      // get copy of shape array
            origin = new int[la.getRank()]; // start with all zeroes
            arrLat = la.read(origin, shape);
            lat_obs = new double[ny_obs];
            index = arrLat.getIndex();
            for (int j = 0; j < ny_obs; j++) {
                index.set(j);
                lat_obs[j] = arrLat.getDouble(index);
            }


        //            lon_obs = (double[]) nc_OBS.findVariable("lon").read(new int[]{0}, new int[]{nx_obs}).reduce().copyToNDJavaArray();
//            lat_obs = (double[]) nc_OBS.findVariable("lat").read(new int[]{0}, new int[]{ny_obs}).reduce().copyToNDJavaArray();


        } catch (IOException e) {
            throw new IOException("Problem reading Obs netcdf file ");
        } catch (InvalidRangeException e) {
            throw new IOException("Problem extracting fields at location " + nc_OBS.getLocation().toString() + " : " +
                    e.getMessage());
        }

        // Recherche de l'indice de début et de fin des années recherchées :
/*
        for (int i = 0; i < liste_year_obs_temp.length; i++) {
        if (liste_year_obs_temp[i] == year_debut) {
        i_deb = i;
        } else if (liste_year_obs_temp[i] == year_fin) {
        i_fin = i;
        }
        }
         */
        for (int ii = 0; ii < lat_obs.length; ii++) {
            int lat_temp = Math.round(Math.round(lat_obs[ii]));
            if (lat_temp == latlim_nord) {
                ilat_deb = ii;
            } else if (lat_temp == latlim_sud) {
                ilat_fin = ii;
            }
        }

        int[] origin_4D = new int[]{0, 0, 0, 0};
        int[] shape_4D = new int[]{NB_YEAR, 12, ny_obs, nx_obs};

        try {


            eggs_obs = (float[][][][]) nc_OBS.findVariable("HueAnch").
                    read(origin_4D, shape_4D).reduce().
                    copyToNDJavaArray();


        } catch (IOException e) {
            throw new IOException("Problem reading Obs 2 output file ");
        } catch (InvalidRangeException e) {
            throw new IOException("Problem extracting fields at location " + nc_OBS.getLocation().toString() + " : " +
                    e.getMessage());
        }
    }

        static void lecture_CLM() throws IOException {
// (la grille est la meme que pour les obs; on ne la relit pas (lon, lat)

        int latlim_nord = -5;
        int latlim_sud = -15;



        for (int ii = 0; ii < lat_obs.length; ii++) {
            int lat_temp = Math.round(Math.round(lat_obs[ii]));
            if (lat_temp == latlim_nord) {
                ilat_deb = ii;
            } else if (lat_temp == latlim_sud) {
                ilat_fin = ii;
            }
        }

        int[] origin_3D = new int[]{0, 0, 0};
        int[] shape_3D = new int[]{12, ny_obs, nx_obs};

        try {
            eggs_CLM = (float[][][]) nc_CLM.findVariable("Hue-anch").
                    read(origin_3D, shape_3D).reduce().
                    copyToNDJavaArray();

        } catch (IOException e) {
            throw new IOException("Problem reading CLM 2 output file ");
        } catch (InvalidRangeException e) {
            throw new IOException("Problem extracting fields at location " + nc_CLM.getLocation().toString() + " : " +
                    e.getMessage());
        }
    }


// -----------------------------------------------------------------------------
    static void get_zonesponte(int year, int month) {
        // recherche les mois et positions (lat lon) on la ponte a été observée, et le stockt dans un tableau {mois lon lat}
        int i_year = 0;
        int i_month = 0;
        i_month = month - 1;
        int i = 0;
        int j = 0;

        lon_ponte = null;
        lat_ponte = null;
        n_ponte = 0;

        for (i = 0; i < liste_year_obs.length; i++) {
            if (liste_year_obs[i] == year) {
                i_year = i;
                break;
            }
        }
//            System.out.println("liste_year_obs[i] =" + liste_year_obs[i] + " ; i_year =" + i_year);
//            System.out.println("month =" + month + " ; i_month = " + i_month);


// 1er passage pour déterminer le nombre de pos ou la ponte a été vue (nb oeuf>0):
        for (i = ilat_deb; i < ilat_fin; i++) {
            for (j = 0; j < nx_obs; j++) {
                if (eggs_obs[i_year][i_month][i][j] > 0 && eggs_obs[i_year][i_month][i][j] != 32767.0f) {
                    //System.out.println("lon = " + lon_obs[j] + " ; lat = " + lat_obs[i] + " ; eggs = " +eggs_obs[i_year][i_month][i][j]);
                    n_ponte++;
                }
            }
        }

        // Declaration de la taille des tableaux :
        lon_ponte = new double[n_ponte];
        lat_ponte = new double[n_ponte];
        densite_ponte = new double[n_ponte];
//        System.out.println("n_ponte =" + Get_SpawnFile.n_ponte);


//        2eme Passage pour remplir les tableaux :
// ATTENTION :      nc{'HueAnch'}.FillValue_ = nclong(32767);
//

        n_ponte = 0;
        for (i = ilat_deb; i < ilat_fin; i++) {
            for (j = 0; j < nx_obs; j++) {
                if (eggs_obs[i_year][i_month][i][j] > 0 && eggs_obs[i_year][i_month][i][j] != 32767.0f) {
                    lon_ponte[n_ponte] = lon_obs[j];
                    lat_ponte[n_ponte] = lat_obs[i];
                    densite_ponte[n_ponte]= eggs_obs[i_year][i_month][i][j];
                    n_ponte++;
                }
            }
        }
    }


    // ----> LIRE ZONE PONTE CLIM :
    static void get_zonesponte_CLM(int month) {
        // recherche les mois et positions (lat lon) on la ponte a été observée, et le stockt dans un tableau {lon lat}
        int i_year = 0;
        int i_month = 0;
        i_month = month - 1;
        int i = 0;
        int j = 0;

        lon_ponte = null;
        lat_ponte = null;
        n_ponte = 0;


// 1er passage pour déterminer le nombre de pos ou la ponte a été vue (nb oeuf>0):
        for (i = ilat_deb; i < ilat_fin; i++) {
            for (j = 0; j < nx_obs; j++) {
                if (eggs_CLM[i_month][i][j] > 0 && eggs_CLM[i_month][i][j] != 32767.0f) {
                    //System.out.println("lon = " + lon_obs[j] + " ; lat = " + lat_obs[i] + " ; eggs = " +eggs_obs[i_year][i_month][i][j]);
                    n_ponte++;
                }
            }
        }
        // Declaration de la taille des tableaux :
        lon_ponte = new double[n_ponte];
        lat_ponte = new double[n_ponte];
        densite_ponte = new double[n_ponte];
//        System.out.println("n_ponte =" + Get_SpawnFile.n_ponte);

//        2eme Passage pour remplir les tableaux :
// ATTENTION :      nc{'HueAnch'}.FillValue_ = nclong(32767);
//

        n_ponte = 0;
        for (i = ilat_deb; i < ilat_fin; i++) {
            for (j = 0; j < nx_obs; j++) {
                if (eggs_CLM[i_month][i][j] > 0 && eggs_CLM[i_month][i][j] != 32767.0f) {
                    lon_ponte[n_ponte] = lon_obs[j];
                    lat_ponte[n_ponte] = lat_obs[i];
                    densite_ponte[n_ponte]= eggs_CLM[i_month][i][j];
                    n_ponte++;
                }
            }
        }
    }

//

    public static void get_zonesponte3(int year, int month) {
        // Cherche les zones de ponte observée depui 1 mois avant jusqu'a 1mois après (sur 3 mois)
        int n_ponte1, n_ponte2, n_ponte3, j;
        double[] lon_ponte1, lon_ponte2, lon_ponte3;
        double[] lat_ponte1, lat_ponte2, lat_ponte3;
        double[] densite_ponte1, densite_ponte2, densite_ponte3;
        int year1, year3, month1, month3;
        year1 = year;
        year3 = year;
        month1 = month - 1;
        month3 = month + 1;

//      System.out.println("year = " + year + "   month = " +  month);
        get_zonesponte(year, month);
        lon_ponte2 = lon_ponte;
        lat_ponte2 = lat_ponte;
        densite_ponte2 = densite_ponte;
        n_ponte2 = n_ponte;

        if (n_ponte > 0) {
            if (month < 2) {
                year1 = year - 1;
                month1 = 12;
            }
            if (month > 11) {
                year3 = year + 1;
                month3 = 1;
            }

//      System.out.println("year1 = " + year1 + "   month1 = " +  month1);
            get_zonesponte(year1, month1);
            lon_ponte1 = lon_ponte;
            lat_ponte1 = lat_ponte;
            densite_ponte1 = densite_ponte;
            n_ponte1 = n_ponte;

//      System.out.println("year3 = " + year3 + "   month3 = " +  month3);
            get_zonesponte(year3, month3);
            lon_ponte3 = lon_ponte;
            lat_ponte3 = lat_ponte;
            densite_ponte3 = densite_ponte;
            n_ponte3 = n_ponte;

            n_ponte_3mois = n_ponte1 + n_ponte2 + n_ponte3;
            lon_ponte_3mois = new double[n_ponte_3mois];
            lat_ponte_3mois = new double[n_ponte_3mois];
            densite_ponte_3mois = new double[n_ponte_3mois];

            for (int i = 0; i < n_ponte1; i++) {
                lon_ponte_3mois[i] = lon_ponte1[i];
                lat_ponte_3mois[i] = lat_ponte1[i];
                densite_ponte_3mois[i] = densite_ponte1[i];
            }
            j = 0;
            for (int i = n_ponte1; i < n_ponte1 + n_ponte2; i++) {
                lon_ponte_3mois[i] = lon_ponte2[j];
                lat_ponte_3mois[i] = lat_ponte2[j];
                densite_ponte_3mois[i] = densite_ponte2[j];
                j++;
            }

            j = 0;
            for (int i = n_ponte1 + n_ponte2; i < n_ponte_3mois; i++) {
                lon_ponte_3mois[i] = lon_ponte3[j];
                lat_ponte_3mois[i] = lat_ponte3[j];
                densite_ponte_3mois[i] = densite_ponte3[j];
                j++;
            }
        } else {
            n_ponte_3mois = 0;
        }
    }

// POUR BOUCHER LES TROU AVEC CLIMATO :
    public static void get_zonesponte3_completeCLM(int year, int month) {
        // Cherche les zones de ponte observée depui 1 mois avant jusqu'a 1mois après (sur 3 mois)
        // SI PAS DE PONTE OBS CE MOIS CI, UTILISE CLIMATO
        int n_ponte1, n_ponte2, n_ponte3, j;
        double[] lon_ponte1, lon_ponte2, lon_ponte3;
        double[] lat_ponte1, lat_ponte2, lat_ponte3;
        double[] densite_ponte1, densite_ponte2, densite_ponte3;

        int year1, year3, month1, month3;
        year1 = year;
        year3 = year;
        month1 = month - 1;
        month3 = month + 1;

//      System.out.println("year = " + year + "   month = " +  month);
        get_zonesponte(year, month);
        lon_ponte2 = lon_ponte;
        lat_ponte2 = lat_ponte;
        densite_ponte2 = densite_ponte;
        n_ponte2 = n_ponte;

        if (n_ponte > 0) { // CAS OU IL Y A DES OBS DE PONTE CE MOIS CI...
            if (month < 2) {
                year1 = year - 1;
                month1 = 12;
            }
            if (month > 11) {
                year3 = year + 1;
                month3 = 1;
            }

//      System.out.println("year1 = " + year1 + "   month1 = " +  month1);
            get_zonesponte(year1, month1);
            lon_ponte1 = lon_ponte;
            lat_ponte1 = lat_ponte;
            densite_ponte1 = densite_ponte;
            n_ponte1 = n_ponte;

//      System.out.println("year3 = " + year3 + "   month3 = " +  month3);
            get_zonesponte(year3, month3);
            lon_ponte3 = lon_ponte;
            lat_ponte3 = lat_ponte;
            densite_ponte3 = densite_ponte;
            n_ponte3 = n_ponte;

            n_ponte_3mois = n_ponte1 + n_ponte2 + n_ponte3;
            lon_ponte_3mois = new double[n_ponte_3mois];
            lat_ponte_3mois = new double[n_ponte_3mois];
            densite_ponte_3mois = new double[n_ponte_3mois];

            for (int i = 0; i < n_ponte1; i++) {
                lon_ponte_3mois[i] = lon_ponte1[i];
                lat_ponte_3mois[i] = lat_ponte1[i];
                densite_ponte_3mois[i] = densite_ponte1[i];

            }
            j = 0;
            for (int i = n_ponte1; i < n_ponte1 + n_ponte2; i++) {
                lon_ponte_3mois[i] = lon_ponte2[j];
                lat_ponte_3mois[i] = lat_ponte2[j];
                densite_ponte_3mois[i] = densite_ponte2[j];
                j++;
            }

            j = 0;
            for (int i = n_ponte1 + n_ponte2; i < n_ponte_3mois; i++) {
                lon_ponte_3mois[i] = lon_ponte3[j];
                lat_ponte_3mois[i] = lat_ponte3[j];
                densite_ponte_3mois[i] = densite_ponte3[j];
                j++;
            }
        } else {
            // CAS OU IL N'Y A PAS DE'OBS DE PONTE CE MOIS CI : COMPLETER AVEC CLIMATO
            get_zonesponte_CLM(month);
            lon_ponte_3mois = lon_ponte;
            lat_ponte_3mois = lat_ponte;
            densite_ponte_3mois = densite_ponte;
            n_ponte_3mois = n_ponte;
        }
    }
}



