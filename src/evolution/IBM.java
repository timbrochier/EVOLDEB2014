package evolution;

public class IBM {

    public static String config_file;
//    public IBM() {
//    }
// -----------------------------------------------------------------------------
    //Main method

    public static void main(String[] args) {
config_file=args[0];
        //Construct the application
        Simulation simulation = new Simulation();

        simulation.initSerial();
        new IBM();
    }
}
      

