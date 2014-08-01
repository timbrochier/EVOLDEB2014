package util;

/**
 *
 * @author timbrochier
 *
 * EN COUR - PEUT ETRE PAS NECESSAIRE (ON SE MET EN "STAND STILL")
 */
public class Evite_cote {
        static boolean cote_devant = false;
        boolean cote_sur_chemin =false;
int i;
float route;
double [] deplacement_cor;


public double[] evitement_cote(double [] position, double [] deplacement){
// ON TOURNE SI YA LA COTE DEVANT :

    //1) y-a-t-il une cote devant dans la direction X?
 //  cote_devant =  check_cote_devant(position, deplacement[0]);
    //2) y-a-t-il une cote devant dans la direction Y?
 //  cote_devant =  check_cote_devant(position, deplacement[1]);

return deplacement_cor;
}
/*
boolean check_cote_devant (double [] position, double deplacement){
boolean cote_dev = false;

        // 1 - Ya-t-il une cote devant, dans la direction X ?
        for ( i=1; i<=nbgridcell_par_dt_adultes; i++){
        route = (float) (i + nbgridcell_par_dt_adultes - Math.floor(nbgridcell_par_dt_adultes));
            // System.out.println("nbgridcell_par_dt_adultes = " + nbgridcell_par_dt_adultes + " ; route = " + route);
        if (dir_banc == 0) {
            cote_devant = !Dataset_EVOL.isInWater(x, y + route) ;
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 1) {
            cote_devant = !Dataset_EVOL.isInWater(x + route, y);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 2) {
            cote_devant = !Dataset_EVOL.isInWater(x, y - route);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 3) {
            cote_devant = !Dataset_EVOL.isInWater(x - route, y);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (i<=nbgridcell_par_dt_adultes && cote_devant){
            cote_sur_chemin=true;}
        }
            cote_devant = (cote_sur_chemin || cote_devant);


  return cote_dev;
}






        // 2 - si il y a une côte devant dans la direction dir_banc, on essaye dir_banc +1, etc.. :
        s = 0;
       while (cote_devant) {
        cote_sur_chemin= false;
            dir_banc = dir_banc + 1;
            if (dir_banc > 3) {
                dir_banc = 0;
            }

            for ( i=1; i<=nbgridcell_par_dt_adultes; i++){

route = (float) (i + nbgridcell_par_dt_adultes - Math.floor(nbgridcell_par_dt_adultes));

        if (dir_banc == 0) {
            cote_devant = !Dataset_EVOL.isInWater(x, y + route);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 1) {
            cote_devant = !Dataset_EVOL.isInWater(x + route, y);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 2) {
            cote_devant = !Dataset_EVOL.isInWater(x, y - route);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (dir_banc == 3) {
            cote_devant = !Dataset_EVOL.isInWater(x - route, y);
//             System.out.println("dir_banc = " + dir_banc + " ; !Dataset_EVOL.isInWater(x, y + i) = " + !Dataset_EVOL.isInWater(x, y + i) + " ; cote_devant = " + cote_devant + " ; i = " +i);
        }
        if (i<=nbgridcell_par_dt_adultes && cote_devant){
            cote_sur_chemin=true;}
        }
            cote_devant = (cote_sur_chemin || cote_devant);

            s ++;
            if (s>5){
      //         System.out.println(" PB : x =  " + x + " ; y = " + y + "  ;  Dataset_EVOL.isInWater(x, y) = " + Dataset_EVOL.isInWater(x, y) + " ; km_par_dt_adultes = " + km_par_dt_adultes);
               System.out.println(" PB : x_debut =  " + x_debut + " ; y = " + y_debut + " ; Dataset_EVOL.isInWater(x_debut, y_debut) = " + Dataset_EVOL.isInWater(x_debut, y_debut) +  " ; nbgridcell_par_jour = " + nbgridcell_par_dt_adultes + " ; dir_banc = " + dir_banc + "  ; cote_sur_chemin = " + cote_sur_chemin + " ; cote_devant = " + cote_devant );
            }
        }

//            System.out.println("dir_banc fin = " + dir_banc);




        Boufe_avant = bouffe;

        if (dir_banc == 0) {
            if (Dataset_EVOL.isInWater(x + uvw[0], y + uvw[1] + nbgridcell_par_dt_adultes)){
                x = x + uvw[0];
                y = y + uvw[1] + nbgridcell_par_dt_adultes;
            }
        } else if (dir_banc == 1) {
             if (Dataset_EVOL.isInWater(x + uvw[0]+ nbgridcell_par_dt_adultes, y + uvw[1] )){
            x = x + uvw[0] + nbgridcell_par_dt_adultes;
            y = y + uvw[1];
            }
        } else if (dir_banc == 2) {
             if (Dataset_EVOL.isInWater(x + uvw[0] , y + uvw[1]- nbgridcell_par_dt_adultes )){
                x = x + uvw[0];
                y = y + uvw[1] - nbgridcell_par_dt_adultes;
                }
        } else if (dir_banc == 3) {
             if (Dataset_EVOL.isInWater(x + uvw[0] - nbgridcell_par_dt_adultes, y + uvw[1])){
            x = x + uvw[0] - nbgridcell_par_dt_adultes;
            y = y + uvw[1];
        }
       }


 *
 */

}
