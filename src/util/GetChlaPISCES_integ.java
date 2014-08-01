package util;


/** import java.io */
import java.io.IOException;
import ucar.ma2.Index;
import java.io.File;


/** import netcdf */
import ucar.ma2.Array;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.dataset.NetcdfDataset;
import evolution.Dataset_EVOL;

public class GetChlaPISCES_integ {


///////////////////////////////
// Declaration of the variables
///////////////////////////////
    /**
     * Grid dimension
     */
    static int nx,  ny;
    /**
 /**
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
     * Chla at current time
     */
    static float[][] Chla_tp0;
    /**
     * /**
     * ... at time t + dt
     */
    static float[][] Chla_tp1;
    /**
     * Time step [days] between two records in NetCDF dataset
     */
    public static int dt_HyMo;
    /**
    /**
     * Name of the Dimension in NetCDF file
     */
    static String strXiDim,  strEtaDim,  strTimeDim;
    /**
     * Name of the Variable in NetCDF file
     */
    static String strChla;
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

        System.out.println("PISCES integ SetUp... ");
        strTimeDim = "TIME";
        strChla = "CHLA_TOT_Z_INTEG";
        strXiDim = "XI_RHO";
        strEtaDim = "ETA_RHO";

        String fileName = get_filename(1992);
        open(fileName);

        getDimNC();
//        readConstantField();

    }

    static String get_filename(int year) {
        // Structure du nom des fichier netcdf roms :
       String fileName = Dataset_EVOL.directory_Suplementary_data + "CHLA_TOT_z_integ50m_ROMSPISCES_Y" + year + ".nc";
        return fileName;
    }

    /**
     * Reads the dimensions of the NetCDF dataset
     * @throws an IOException if an error occurs while reading the dimensions.
     */
    private static void getDimNC() throws IOException {

        try {
            nx = ncIn.findDimension(strXiDim).getLength();
            ny = ncIn.findDimension(strEtaDim).getLength();
            //nz = ncIn.findDimension(strZDim).getLength();
        } catch (NullPointerException e) {
            throw new IOException("Problem reading dimensions from dataset " + ncIn.getLocation() + " : " + e.getMessage());
        }
        System.out.println("nx = " + nx);
        System.out.println("ny = " + ny);

        ipo = jpo = 0;
    }

    public static void init(int year) throws IOException {
        System.out.println("initialisation... ");
        //long t0 = Simulation.get_t0();
        // ici on ne repere plus les fichiers par leur "t0" mais par leur nom directement:
        // open(getFile(t0));
        String fileName = get_filename(year);
        System.out.println("fileName : " + fileName);

        open(fileName);
System.out.println("strTimeDim : " + strTimeDim);
//        nbTimeRecords = ncIn.findDimension(strTimeDim).getLength();

        int rank;
        setAllFieldsTp1AtTime(rank = 0);
    }



    /**
     * Reads time dependant variables in NetCDF dataset at specified rank.
     * @param rank an int, the rank of the time dimension in the NetCDF dataset.
     * @throws an IOException if an error occurs while reading the variables.
     */
    static void setAllFieldsTp1AtTime(int rank) throws IOException {

        // int[] origin = new int[]{rank, 0, jpo, ipo};
        //double time_tp0 = time_tp1;

        try {

            Chla_tp1 = (float[][]) ncIn.findVariable(strChla).read(
                    new int[]{rank, 0, 0},
                    new int[]{1, ny, nx}).reduce().copyToNDJavaArray();

    //System.out.println("Chla_tp1[138][185] = " + Chla_tp1[185][138]);


        } catch (IOException e) {
            throw new IOException("Problem 1 extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (InvalidRangeException e) {
            throw new IOException("Problem 2 extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        } catch (NullPointerException e) {
            throw new IOException("Problem 3 extracting fields at location " + ncIn.getLocation().toString() + " : " +
                    e.getMessage());
        }
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

        try {if (ncIn == null ||
                    (new File(ncIn.getLocation()).compareTo(new File(filename)) !=
                    0)) {

//                MainFrame.getStatusBar().setMessage(Resources.MSG_OPEN +
//                        filename);
                System.out.println("On ouvre le fichier : " + filename);
                ncIn = NetcdfDataset.openFile(filename, null);
                nbTimeRecords =
                        ncIn.findDimension(strTimeDim).getLength();
                        System.out.println("nbTimeRecords CHLA_PISCES_INTEG: " + nbTimeRecords );

//            System.out.print("POINT 2 Open dataset " + filename + "\n");
}
        } catch (IOException e) {
            throw new IOException("Problem opening dataset " + filename + " - " + e.getMessage());
        } catch (NullPointerException e) {
            throw new IOException("Problem reading " + strTimeDim + " dimension at location " + filename +
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
    public static void loadChla_PISCES_integ(int jour, int year) throws IOException {

        // OUverture du fichier nc (se fait seulement s'il n'est pas encore ouvert en principe)
        String fileName = get_filename(year);
        open(fileName);

        // ICI ON CONSIDERE QUE :
        //  - LE NOMBRE D'ENREGISTREMENT ROMS PAR MOIS EST CONSTANT POUR LES MOIS DE 1 A 11 (nbTimeRecords);
        //  - LE DERNIER MOIS PEUT EVENTUELEMENT AVOIR UN NOMBRE DIFFERENT D'ENREGISTREMENT (nbTimeRecords)
        int month = 1 + ((int) jour / 30) % 12; // ON CONSIDERE QUE LES MOIS ONT 30 JOURS...
        int jourdumois = (int) jour % 30;
        if (jour == 360) {
            jourdumois = 30;
            month = 12;
        } // pour le cas ou on a 365 jours...

        int rank = month-1;
        // 1) COPIER LES DONNÉE PRécédantes "Tp1" dans "Tp0"
        //setAllFieldsAtTime();
        setAllFieldsTp1AtTime(rank);
    }

public static double getCHLA_PISCES_integ(double xRho, double yRho) throws
        ArrayIndexOutOfBoundsException {
        double chla;
        //-----------------------------------------------------------
        // Interpolate the Chla field
        // in the computational grid.
        int i = (int) Math.round(Math.round(xRho));
        int j = (int) Math.round(Math.round(yRho));

        chla = Chla_tp1[j][i];
        return chla;
    }

public static double getCHLA_PISCES_integ(int i, int j) throws
            ArrayIndexOutOfBoundsException {
        double chla;
        chla = Chla_tp1[j][i];
        return chla;
    }

    //---------- End of class
}