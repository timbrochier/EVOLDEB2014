package util;

 /**
    * * Classe : Configuration
    * *
    * * Description :
    * * La classe Configuration va obtenir et mettre à jour les infos de
    * * configuration à partir des fichiers de configuration
    * * stockés sous [rep_courant]/config
    * * Le fichier doit existé ou une exception sera levée.
    * *
    * * date : 8/05/2004
    * * @author : PeeX Team
    * * @version : 1.0
    * *
    * */
    // Les import

     import java.io.FileInputStream;
     import java.io.FileOutputStream;
     import java.util.Properties;
     import java.lang.Exception;
     import java.io.IOException;

     public class Lire_configuration{
     /**
    * * Methode : getConfig
    * *
    * * Description :
    * * La méthode getConfig va retourner l'information de configuration
    * * désirée à partir d'un fichier de configuration
    * *
    * * date : 8/05/2004
    * * @param : String fichier : Le nom du fichier de configuration
    * * @param : String key : La clé dont on veut obtenir la valeur
    * *
    * * @return : String représentant la valeur de l'info
    * *
    * */

    static Properties config;

         public static void Setup(String fichier)throws IOException {
//     System.out.println("Fichier config = " + fichier);
//     System.out.println("Key = " + key);
     // On fait pointer notre Properties sur ke fichier
     FileInputStream fis = new FileInputStream(fichier);
     config = new Properties();
     config.load(fis);

     fis.close();
     // C'est important de mettre à null, le garbage collector
     // passe plus vite !
///     config = null;
     fis = null;

         }


     public static String getConfig_string(String key) throws Exception
     {

     String tmp = config.getProperty(key);
     return tmp;

     }

     public static boolean getConfig_boolean(String key) throws Exception
     {
        boolean tmp = Boolean.parseBoolean(config.getProperty(key));
     return tmp;

     }

          public static int getConfig_int(String key) throws Exception
     {
     int tmp = Integer.parseInt(config.getProperty(key));
     return tmp;

     }

          public static double getConfig_double(String key) throws Exception
     {
     double tmp = Double.parseDouble(config.getProperty(key));
     return tmp;

     }
     /**
    * * Methode : setConfig
    * *
    * * Description :
    * * La méthode setConfig va mettre à jour/ inserer l'information de configuration
    * * désirée à partir dans un fichier de configuration
    * *
    * * date : 8/05/2004
    * * @param : String fichier : Le nom du fichier de configuration
    * * @param : String key : La clé dont on veut obtenir la valeur
    * * @param : String valeur : La valeur associée à la clé
    * *
    * * @return : String représentant la valeur de l'info
    * *
    * */
     public void setConfig(String fichier, String key, String valeur) throws Exception
     {
     // La petite feinte : Il faur recharger entièrement le fichier
     // et le réecrire.

     //On construit l'adresse du fichier
     String leFichier = System.getProperty("user.dir") + "/config/" + fichier;

     // On fait pointer notre Properties sur le fichier
     Properties config = new Properties();
     FileInputStream fis = new FileInputStream(leFichier);
     config.load (fis);
     fis.close();
     FileOutputStream fos = new FileOutputStream(leFichier);

     config.setProperty(key,valeur);

     config.store (fos,"Dernière mise a jour :");
     // C'est important de mettre à null, le garbage collector
     // passe plus vite !
     fos.close();
     leFichier = null;
     fos = null;
     fis = null;
     config = null;
     }
     }
