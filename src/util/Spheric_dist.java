package util;

import evolution.Simulation;

/**
 *
 * @author timbrochier
 */
public class Spheric_dist {

    double rayon_terre = 6371.01; // en Km
    static double delta_lon, delta_lat, dist;
    
            /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
            /*::                                                                         :*/
            /*::  This routine calculates the distance between two points (given the     :*/
            /*::  latitude/longitude of those points). It is being used to calculate     :*/
            /*::  the distance between two ZIP Codes or Postal Codes using our           :*/
            /*::  ZIPCodeWorld(TM) and PostalCodeWorld(TM) products.                     :*/
            /*::                                                                         :*/
            /*::  Definitions:                                                           :*/
            /*::    South latitudes are negative, east longitudes are positive           :*/
            /*::                                                                         :*/
            /*::  Passed to function:                                                    :*/
            /*::    lat1, lon1 = Latitude and Longitude of point 1 (in decimal degrees)  :*/
            /*::    lat2, lon2 = Latitude and Longitude of point 2 (in decimal degrees)  :*/
            /*::    unit = the unit you desire for results                               :*/
            /*::           where: 'M' is statute miles                                   :*/
            /*::                  'K' is kilometers (default)                            :*/
            /*::                  'N' is nautical miles                                  :*/
            /*::  United States ZIP Code/ Canadian Postal Code databases with latitude & :*/
            /*::  longitude are available at http://www.zipcodeworld.com                 :*/
            /*::                                                                         :*/
            /*::  For enquiries, please contact sales@zipcodeworld.com                   :*/
            /*::                                                                         :*/
            /*::  Hexa Software Development Center © All Rights Reserved 2003            :*/
            /*::                                                                         :*/
            /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
            // Modif/corrections : Timothée
    public

    static double distance(double lat1, double lon1, double lat2, double lon2) {
        delta_lon = lon1 - lon2;
        dist = Math.sin(deg2rad(lat1)) * Math.sin(deg2rad(lat2)) +
                Math.cos(deg2rad(lat1)) * Math.cos(deg2rad(lat2)) *
                Math.cos(deg2rad(delta_lon));
        dist = Math.acos(dist); //<-- c'est la distance angulaire
        //dist = rad2deg(dist);
        //dist = dist * 60 * 1.1515;
        if (Simulation.dist_unit.equals("Km")) {
            dist = dist * 6371.01;
        } else if (Simulation.dist_unit.equals("Nm")) {
            dist = dist * 3440.07;
        }
        return (dist);
    }
// Méthode + précise pour les tres faibles distances (<1km) :
    public static double distance_haversine(double lat1, double lon1, double lat2, double lon2) {
        delta_lon = lon1 - lon2;
        delta_lat = lat1 - lat2;
        dist = Math.sin(deg2rad(delta_lat/2))*Math.sin(deg2rad(delta_lat/2)) +
                Math.cos(deg2rad(lat1)) * Math.cos(deg2rad(lat2)) *
                Math.sin(deg2rad(delta_lon/2))*Math.sin(deg2rad(delta_lon/2));
        dist = 2*Math.asin(Math.sqrt(dist)); //<-- c'est la distance angulaire
        //dist = rad2deg(dist);
        //dist = dist * 60 * 1.1515;
        if (Simulation.dist_unit.equals("Km")) {
            dist = dist * 6371.01;
        } else if (Simulation.dist_unit.equals("Nm")) {
            dist = dist * 3440.07;
        }
        return (dist);
    }


   // Distance eulérienne (NON SPHERIQUE)
    public static double distance_euler(double lat1, double lon1, double lat2, double lon2) {
        delta_lon = lon1 - lon2;
        delta_lat = lat1 - lat2;
        dist = delta_lon*delta_lon + delta_lat*delta_lat;
        dist = Math.sqrt(dist); //<-- c'est la distance

        if (Simulation.dist_unit.equals("Km")) {
            dist = dist * 110; // (car 110 Km à l'équateur
        } else if (Simulation.dist_unit.equals("Nm")) {
            dist = 0000;
        }
        return (dist);
    }



    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    /*::  This function converts decimal degrees to radians             :*/
    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    private static double deg2rad(double deg) {
        return (deg * 2 * Math.PI / 360.0);
    }

    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    /*::  This function converts radians to decimal degrees             :*/
    /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
    private double rad2deg(double rad) {
        return (rad * 360 / (2 * Math.PI));
    }
}