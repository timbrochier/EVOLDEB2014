package util;

import evolution.Simulation;

/**
 * <p>This class attempts to simulate the growth, in length, of ichthyoplankton.
 * By default it simulates the growth as a function of water temperature.
 * It also gives the possibility to use outputs of physical-biogeochemical
 * coupled models like ROMS-BIO (Koné et al. 2005). This allows
 * considering growth of feeding larvae as functions of both temperature and
 * prey availability.</p>
 * The predefined early life stages are egg, yolk-sac larva and feeding larva.
 * Particles change stages at predefined lengh thresholds defined above as
 * constants. All the parameters dealing with the growth model are defined
 * directly in the source code and they are not accessible from the
 * configuration editor. Growth models are specific for each species, and the
 * generic purpose of this application does not allow to provide general
 * parametrization. Therefore, this class must be seen as an example of what
 * can be done in terms of growth modeling. Each biologist will then have to
 * adapt the parameters and the growth function.
 * <p>
 * Reference:
 * Koné, V. (2006, in french). Modélisation de la production primaire et
 * secondaire de l'écosystème du Benguela sud.
 * Influence des conditions trophiques sur le recrutement des
 * larves d'anchois. Ph.D. thesis, Université Pierre & Marie Curie (Paris VI).
 * </p>

 * @author P.Verley
 */
public class GrowthModel {

///////////////////////////////
// Declaration of the constants
///////////////////////////////
    /**
     * Characterized the egg stage. In this model, an egg is particle with a
     * length smaller to the hatch-length defined below.
     */
    final public static int EGG = 0;
    /**
     * Characterized the "yolk sac" stage. A yolk sac larva has a length ranging
     * from the hatch-length and the yolk-to-feeding-length defined below.
     */
    final public static int YOLK_SAC_LARVA = 1;
    /**
     * Characterized the "feeding larva" stage. A feeding larva has a length
     * bigger than the yolk-to-feeding-length defined below.
     */
    final public static int FEEDING_LARVA = 2;
    /**
     * The growth function assumed the sea water temperature must not be
     * be colder than this threshold. Temperature set in Celsius degree.
     */
    final private static double TP_THRESHOLD = 10.d; //°C
    /**
     * Initial length [millimeter] for the particles.
     */
    final public static double LENGTH_INIT = 0.025d; // mm
    /**
     * Threshold [millimiter] to distinguish eggs from larvae
     */
    final public static double HATCH_LENGTH = 2.8d; //mm
    /**
     * Threshold [millimeter] between Yolk-Sac Larvae and Feeding Larvae
     */
    final private static double YOLK_TO_FEEDING_LENGTH = 4.5d; //mm
    /**
     * Half saturation constant.
     */
    final private static float KS = 0.5f;
    /**
     * Preferency for large Diatoms [0 ; 1]
     */
    final private static float E21 = 1.f / 4.f;
    /**
     * Preferency for NanoPhyto [0 ; 1]
     */
    final private static float E22 = 1.f / 4.f;
    /**
     * Preferency for MicroZoo [0 ; 1]
     */
    final private static float E23 = 1.f / 4.f;
    /**
     * Preferency for MesoZoo [0 ; 1]
     */
    final private static float E24 = 1.f / 4.f;

///////////////////////////////
// Declaration of the variables
///////////////////////////////
    /**
     * <code>true</code> if dead cold, <code>false</code> otherwise.
     */
    private static boolean deadCold;
    /**
     * The stage of the particle: egg, yolk-sac larva or feeding larva.
     * @see ichthyop.util.Constant for labels characterizing the stages.
     */
    private static int stage;
    /**
     * Lethal see water temperature [celsius] for egg
     */
    static double lethalTpEgg;
    /**
     * Lethal see water temperature [celsius] for larva
     */
    static double lethalTpLarvae;
    /**
     * Model time step expressed in days.
     */
    private static double dt_day;
    /**
     * Flag that indicates whether lethal water temperature should be simulated.
     */
    private static boolean FLAG_LETHAL_TP;
    /**
     * Estimation of the prey availability.
     */
    private static double food;
    /**
     * Growth limitation factor function of prey availability.
     */
    private static double foodLimFactor;

////////////////////////////
// Definition of the methods
////////////////////////////
    /**
     * Initializes the growth model. Gets the values of required parameters.
     */
    public static void init() {

        FLAG_LETHAL_TP = true; // (On considere toujours une temperature letale)
        if (FLAG_LETHAL_TP) {
            // Pour le moment on met la même temp letale pour oeuf et larves
            lethalTpEgg = Simulation.temp_lethal_min;
            lethalTpLarvae = Simulation.temp_lethal_min;
        }
        deadCold = false;
        // tim : dt en heures
        dt_day = (double) Simulation.dt_advec / (24*60); // (24heures/jour)

    }

    /**
     * Makes the particle grow, as a linear function of water temperature.
     *
     * @param length a double, the current length of the particle.
     * @param temperature a double, the sea water temperature [celsius] at
     * particle location.
     * @return a double, the new length of the particle.
     */
    public static double grow(double length, double temperature) {

        stage = getStage(length);
        deadCold = FLAG_LETHAL_TP
                ? ((stage == EGG) && (temperature < lethalTpEgg)) || ((stage > EGG) && (temperature < lethalTpLarvae))
                : false;
        if (!deadCold) {
            length += (.02d + .03d * Math.max(temperature, TP_THRESHOLD)) *
                    dt_day;
            return length;
        }
        return length;
    }

    /**
     * Makes the particle grow, as both function of water temperature and prey
     * availability.
     *
     * @param length a double, the current length of the particle.
     * @param temperature a double, the sea water temperature [celsius] at
     * particle location.
     * @param lPhyto a double, concentration [mMol/m3] in large phytoplankton
     * at particle location.
     * @param sZoo a double, concentration [mMol/m3] in small zooplankton
     * at particle location.
     * @param lZoo a double, concentration [mMol/m3] in large zooplankton
     * at particle location.
     * @return a double, the new length of the particle.
     */
    public static double grow(double length, double temperature, double Diatoms,
            double NanoPhyto, double MicroZoo, double MesoZoo) {

        stage = getStage(length);
        deadCold = FLAG_LETHAL_TP
                ? (((stage == EGG) && (temperature < lethalTpEgg)) || ((stage > EGG) && (temperature < lethalTpLarvae)))
                : false;
        if (!deadCold) {
            switch (stage) {
                case EGG:
                case YOLK_SAC_LARVA:
                    foodLimFactor = 1.f;
                    break;
                case FEEDING_LARVA:
                    // Dans PISCES, le plankton est en uMol/L --> equivalent aux mMol/m3
                    food = E21 * Diatoms + E22 * NanoPhyto + E23 * MicroZoo + E24 * MesoZoo;
                    foodLimFactor = food / (KS + food);
                    //System.out.println("foodLimFactor = " + foodLimFactor);
                    break;
            }
            length += foodLimFactor * (.02d + .03d * Math.max(temperature,
                    TP_THRESHOLD)) * dt_day;
        }
        return length;
    }

    /**
     * Gets the stage of the particle, given the specified length.
     *
     * @param length a double, the length [millimeter] of the particle.
     * @return an int that characterized the stage of the particle (egg,
     * yolk-sac larva or feeding larva).
     * @see ichthyop.util.Constant for details about the labels characterizing
     * the stage of the particle.
     */
    public static int getStage(double length) {

        /** Yolk-Sac Larvae */
        if (length >= HATCH_LENGTH & length < YOLK_TO_FEEDING_LENGTH) {
            return YOLK_SAC_LARVA;
        } /** Feeding Larvae */
        else if (length >= YOLK_TO_FEEDING_LENGTH) {
            return FEEDING_LARVA;
        }
        /** eggs */
        return EGG;
    }

    /**
     * Determines whether the particle is dead cold.
     *
     * @return <code>true</code> if dead cold, <code>false</code> otherwise.
     */
    public static boolean isDeadCold() {
        return deadCold;
    }

    //----------- End of class
}
