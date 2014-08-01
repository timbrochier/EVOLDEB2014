package util;


/** import java.io */
import java.io.IOException;
import ucar.ma2.Index;


/** import netcdf */
import ucar.ma2.Array;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.dataset.NetcdfDataset;
import evolution.Dataset_EVOL;

public class GetOxyclinDepth_BIDOUILLE {


///////////////////////////////
// Declaration of the variables
///////////////////////////////
    /**
     * Grid dimension
     */
    static int nx,  ny;

    // liste des années
    static int[] liste_year;
 /**
     *

     * Origin for grid index
     */
    static int ipo,  jpo;
    /**
     * Number of time records in current NetCDF file
     */
    //static int nbTimeRecords_init; // pour le premier mois
    static int nbTimeRecords; // Pour vérification chaque moi (car change dans certaines sorties ROMS)
    /**
     * The NetCDF dataset
     */
    static NetcdfFile ncIn;
    /**
     * Longitude at rho point.
     */
    static double[][] lonRho;
    /**
     * Latitude at rho point.
     */
    static double[][] latRho;
    /**
     * Bathymetry
     */
    static double[][] hRho;

    /**
     * OxyClin at current time
     */
    static float[][] Oxyclin_tp0;
    /**
     * /**
     * ... at time t + dt
     */
    static float[][] Oxyclin_tp1;
    /**
     * Geographical boundary of the domain
     */
    private static double latMin,  lonMin,  latMax,  lonMax;
    /**
     * Time step [days] between two records in NetCDF dataset
    /**
     * Name of the Dimension in NetCDF file
     */
    static String strXiDim,  strEtaDim;
    /**
     * Name of the Variable in NetCDF file
     */
    static String strOxyclin, strYearDim, strDtDim;
    /**
     * Name of the Variable in NetCDF file
     */
    static String strLon,  strLat;
    /**


///////////////////////////////
// Declaration of the constants
///////////////////////////////
    ////////////////////////////
// Definition of the methods
////////////////////////////
    /**
     * Sets up the {@code Dataset}. The method first sets the appropriate
     * variable names, loads the first NetCDF dataset and extract the time
     * non-dependant information, such as grid dimensions, geographical
     * boundaries, depth at sigma levels.
     * @throws an IOException if an error occurs while setting up the
     * {@code Dataset}
     */
    public static void setUp() throws IOException {

        strYearDim = "year";
        strDtDim = "dt_avg_5j";
        strXiDim = "xi_rho";
        strEtaDim = "eta_rho";
        strLon = "lon_rho";
        strLat = "lat_rho";
        strOxyclin = "Oxy_1ml_BIDOUILLE";


        String fileName = get_filename();
        open(fileName);

        getDimNC();
        readConstantField();
        getDimGeogArea();

    }

    static String get_filename() {
        // Structure du nom des fichier netcdf roms :
       String fileName = Dataset_EVOL.directory_Suplementary_data + Dataset_EVOL.oxyclin_filename;
        return fileName;
    }

    /**
     * Reads the dimensions of the NetCDF dataset
     * @throws an IOException if an error occurs while reading the dimensions.
     */
    private static void getDimNC() throws IOException {
            int NB_YEAR = 0; // number of year with observations
        int[] shape;
        int[] origin;
        Array arrYear;
        Index index;

        try {
            NB_YEAR = ncIn.findDimension("year").getLength();

            nx = ncIn.findDimension(strXiDim).getLength();
            ny = ncIn.findDimension(strEtaDim).getLength();
            //nz = ncIn.findDimension(strZDim).getLength();
               System.out.println("nx = " + nx);
        System.out.println("ny = " + ny);
        //System.out.println("nz = " + nz);

        ipo = jpo = 0;

// GET LIST YEAR :
            ucar.nc2.Variable y = ncIn.findVariable("year");
            shape = y.getShape();      // get copy of shape array
            origin = new int[y.getRank()]; // start with all zeroes
            arrYear = y.read(origin, shape); //.reduce().copyToNDJavaArray();    // read data from file
            liste_year = new int[NB_YEAR];
            index = arrYear.getIndex();
            for (int j = 0; j < NB_YEAR; j++) {
                index.set(j);
                liste_year[j] = arrYear.getInt(index);
            //   System.out.println("liste_year_obs[j] = " + liste_year_obs[j]);
            }

 } catch (NullPointerException e) {
            throw new IOException("Problem reading dimensions from dataset " + ncIn.getLocation() + " : " + e.getMessage());
        } catch (IOException e) {
            throw new IOException("Problem reading Obs netcdf file ");
        } catch (InvalidRangeException e) {
            throw new IOException("Problem extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        }
        }

    /**
     * Determines the geographical boundaries of the domain in longitude,
     * latitude and depth.
     */
    private static void getDimGeogArea() {

        //--------------------------------------
        // Calculate the Physical Space extrema

        lonMin = Double.MAX_VALUE;
        lonMax = -lonMin;
        latMin = Double.MAX_VALUE;
        latMax = -latMin;
        int i = nx;
        int j = 0;

        while (i-- > 0) {
            j = ny;
            while (j-- > 0) {
                if (lonRho[j][i] >= lonMax) {
                    lonMax = lonRho[j][i];
                }
                if (lonRho[j][i] <= lonMin) {
                    lonMin = lonRho[j][i];
                }
                if (latRho[j][i] >= latMax) {
                    latMax = latRho[j][i];
                }
                if (latRho[j][i] <= latMin) {
                    latMin = latRho[j][i];
                }
            }
        }
        //System.out.println("lonmin " + lonMin + " lonmax " + lonMax + " latmin " + latMin + " latmax " + latMax);
        //System.out.println("depth max " + depthMax);

        double double_tmp;
        if (lonMin > lonMax) {
            double_tmp = lonMin;
            lonMin = lonMax;
            lonMax = double_tmp;
        }

        if (latMin > latMax) {
            double_tmp = latMin;
            latMin = latMax;
            latMax = double_tmp;
        }
    }

    /**
     * Initializes the {@code Dataset}. Opens the file holding the first time
     * of the simulation. Checks out the existence of the fields required
     * by the current simulation. Sets all fields at time for the first time
     * step.
     * @throws an IOException if a required field cannot be found in the NetCDF
     * dataset.
     */
    public static void init() throws IOException {
        System.out.println("initialisation... ");
        //long t0 = Simulation.get_t0();
        // ici on ne repere plus les fichiers par leur "t0" mais par leur nom directement:
        // open(getFile(t0));
        String fileName = get_filename();
       // System.out.println("fileName : " + fileName);

        open(fileName);
//        nbTimeRecords = ncIn.findDimension(strTimeDim).getLength();

        setAllFieldsTp1AtTime(0,0);
    }

    /**
     * Opens the NetCDF dataset from the specified location.
     * If the <code>rawPath</code> is an OPeNDAP URL, it directly opens it.
     * If the <code>rawPath</code> is a local path, the application first lists
     * the files of the folder by a call to {@code getInputList} method and then
     * opens the first file of the list.
     *
     * @param rawPath a String, the location of the dataset in URI format.
     * It can be local path or an OPeNDAP URL.
     * @throws an IOException if an erroc occurs when opening the dataset.
     * @see java.net.URI for details about URI syntax
     */

    /**
     * Loads the NetCDF dataset from the specified filename.
     * @param filename a String that can be a local pathname or an OPeNDAP URL.
     * @throws IOException
     */
    private static void open(String filename) throws IOException {

        try {
//                MainFrame.getStatusBar().setMessage(Resources.MSG_OPEN +
//                        filename);
                System.out.println("On ouvre le fichier : " + filename);
                ncIn = NetcdfDataset.openFile(filename, null);


        } catch (IOException e) {
            throw new IOException("Problem opening dataset " + filename + " - " + e.getMessage());
        } catch (NullPointerException e) {
            throw new IOException("Problem reading dimension at location " + filename +
                    " : " + e.getMessage());
        }

    }

    /**
     * Reads time non-dependant fields in NetCDF dataset
     */
    static void readConstantField() throws IOException {

        int[] origin = new int[]{jpo, ipo};
        int[] size = new int[]{ny, nx};
        Array arrLon, arrLat;
        Index index;

        StringBuffer list = new StringBuffer(strLon);
        list.append(", ");
        list.append(strLat);
        list.append(", ");
        try {
            arrLon = ncIn.findVariable(strLon).read(origin, size);
            arrLat =
                    ncIn.findVariable(strLat).read(origin, size);

            if (arrLon.getElementType() == double.class) {
                lonRho = (double[][]) arrLon.copyToNDJavaArray();
                latRho = (double[][]) arrLat.copyToNDJavaArray();
            } else {
                lonRho = new double[ny][nx];
                latRho =
                        new double[ny][nx];
                index =
                        arrLon.getIndex();
                for (int j = 0; j <
                        ny; j++) {
                    for (int i = 0; i <
                            nx; i++) {
                        index.set(j, i);
                        lonRho[j][i] = arrLon.getDouble(index);
                        latRho[j][i] = arrLat.getDouble(index);
                    }
                }
            }

        } catch (IOException e) {
            throw new IOException("Problem (1) reading one of the fields " + list.toString() + " at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (InvalidRangeException e) {
            throw new IOException("Problem (2) reading one of the fields " + list.toString() + " at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (NullPointerException e) {
            throw new IOException("Problem (3) reading one of the fields " + list.toString() + " at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        }

    }


// Timothee...
    public static void load_Oxyclin(int jour, int year) throws IOException {

        int i_year = 0;
        int i_jour = 0;
        // i_jour = jour - 1;   A VERIFIER
        int i = 0;
        int j = 0;

        for (i = 0; i < liste_year.length; i++) {
            if (liste_year[i] == year) {
                i_year = i;
                break;
            }
        }

        // i_year est l'indice de l'annee

        i_jour = jour/Dataset_EVOL.dt_HyMo;

        // i_jour est l'indice du temps pour une année donne

        // 1) COPIER LES DONNÉE PRécédantes "Tp1" dans "Tp0"
        //setAllFieldsAtTime();
        setAllFieldsTp1AtTime(i_year,i_jour);
    }


    /**
     * Reads time dependant variables in NetCDF dataset at specified rank.
     * @param rank an int, the rank of the time dimension in the NetCDF dataset.
     * @throws an IOException if an error occurs while reading the variables.
     */
    static void setAllFieldsTp1AtTime(int i_year, int i_jour) throws IOException {

//        int[] origin = new int[]{i_jour, i_year, jpo, ipo};

        try {

            Oxyclin_tp1 = (float[][]) ncIn.findVariable(strOxyclin).read(
                    new int[]{i_jour,i_year, 0, 0},
                    new int[]{1,1, ny, nx}).reduce().copyToNDJavaArray();


        } catch (IOException e) {
            throw new IOException("Problem extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (InvalidRangeException e) {
            throw new IOException("Problem extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (NullPointerException e) {
            throw new IOException("Problem extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        }
    }


public static double getOxyclin_depth(double xRho, double yRho) throws
            ArrayIndexOutOfBoundsException {
        double Oxy;
        //-----------------------------------------------------------
        // Interpolate the Chla field
        // in the computational grid.
        int i = (int) Math.round(Math.round(xRho));
        int j = (int) Math.round(Math.round(yRho));

        if (Oxyclin_tp1[j][i] == 32767.0){
           Oxy = -99999;
        }else {
          Oxy = -Oxyclin_tp1[j][i];

        }
        return Oxy;
    }

public static double getOxyclin_depth(int i, int j) throws
            ArrayIndexOutOfBoundsException {
        double Oxy;

        if (Oxyclin_tp1[j][i] == 32767.0){
           Oxy = -99999;
        }else {
          Oxy = -Oxyclin_tp1[j][i];

        }
        return Oxy;
    }
    //---------- End of class
}

