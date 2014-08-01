package util;

/** import netcdf data */
/**
 * Both ROMS and MARS use an Arakawa C grid.
 * The Rho point represents a 2D or 3D point within the C grid. It is referenced
 * by two sets of coordinates:
 * <ul>
 * <li>grid coordinates (x, y, z) if 3D, (x, y) if 2D
 * <li>geographical coordinates (longitude, latitude) if 3D,
 * (longitude, latitude) if 2D
 * </ul>
 * The class provides methods to switch from grid coordinates (x, y, z)
 * to geographical coordinates (lon, lat, depth) and reciprocally.
 *
 * @author P.Verley
 */
public class Position_ponte {



///////////////////////////////
// Declaration of the variables
///////////////////////////////


    public final int date_grid;
	public final int longitude_grid;
	public final int latitude_grid;
	public final int profontdeur_grid; // ATTENTION, indice de la liste des prof portentielles, pas sigma


///////////////
// Constructors
///////////////


   public Position_ponte(int date_grid, int longitude_grid, int latitude_grid, int profontdeur_grid) {

        this.date_grid = date_grid;
		this.longitude_grid = longitude_grid;
		this.latitude_grid = latitude_grid;
		this.profontdeur_grid = profontdeur_grid;
	}




////////////////////////////
// Definition of the methods
////////////////////////////

// ...

//////////
// Setters
//////////

//...

//////////
// Getters
//////////

//...

    //----------- End of class
}

