
            // test pas de temps variable
            // (il faut initialiser les old sur la base d'un pas de temps de 1h (arbitraire mais intermediare)
            // (REMPLACER PAR LE TEMP QUE UN POISSON DE 20 CM MET A PARCOURIR 10 KM EN VITESSE DE CROISIERE)
             // 0.2 m/s * 3600s *Vitesse_max_Bls(8)=  5.7600 km par heure
            //--> pour 10km dt = 10km / 5.7600km/h = 1.7361 heure
            //                  = round(1.7361 * 60) minutes = 104 minutes
            // Si dt < 1 minute (ex :Bls = 20)
            // dt = 1minute (au pire ils avancent de moins d'une minute de trop)
/*
            //distance parourue, en unite grille, durant le dernier pas de temps : dist_grd_old
            double dist_old_km = Math.sqrt(Math.pow(uvw[0],2) + Math.pow(uvw[1],2) + Math.pow(uvw[2],2))*Population.mean_grid_size;
            double dist_old_km = util.Spheric_dist.distance_haversine(Dataset_EVOL.lat_min, Dataset_EVOL.lon_min, Dataset_EVOL.lat_max, Dataset_EVOL.lon_max);

            // dernier pas de temps = dt_old;
            double Vt_old = dist_old/dt_adapt_old; //vitesse precedante (en unite grille par minute) (dt_adapt = minutes)

            dt_present = 10km/Vt_old; // Vt_old en [km/h]
            frac = frac_old + 1/dt_present;

            // on enregistre les old :
            frac_old = frac;
            Vt_old = Vt;
            dt_old = dt;
*/
            // fin test


        // DANS CETTE VERSION : pas d'aspect cognitif,
        // le banc avance toujours droit devant lui, sauf si les conditions ne sont plus "bonnes"
        // la il tourne et il avance encore
        // la direction de déplacement est donc une propriété du banc :
        // dir_banc = 0,1,2,3 (Nord, est, sud, ouest)
        // conditions de phyto ici : Boufe_avant

        double[] SST_voisin = new double[4];
        int[] BAT_voisin = new int[4];




// A DEFINIR EN FONCTION DU DEB :
        Body_length = DEB.getLength(); // cm


        // Valeures de vit_banc_Bls : (These P. Brehmer)
            // "vitesse d'exploration" = 2
            // "vitesse max = 8
            // "vitesse moyene = 4

        //--> 1/Température (SST) ici est-elle dans l'optimal temperature range?
        // Si oui on ne bouge pas, sinon on regarde si c'est mieux autour
        float optimal_temp_min, optimal_temp_max;
        optimal_temp_min = 16;
        optimal_temp_max = 22;        //22; // Le Fur et al., 2009 : T°C opt = 21-25

//               System.out.println("Phyto ICI = " + plancton[1]);


// ON TOURNE SI LES CONDITIONS SONT MOINS BONNES QU'AVANT :
        double bouffe_f = bouffe/ (bouffe+ DEB.X_K);

        // if ((temp < optimal_temp_min || temp > optimal_temp_max) ||
//            System.out.println("dir_banc début = " + dir_banc);

        if (bouffe_f<0.5){
            dir_banc = 1; // (Est)
            vit_banc_Bls = 1;
        //System.out.println("CE POISSON RANDOM FONCE VERS L'EST CAR IL A VRAIMENT PLUS RIEN A BECTER (bouffe_f<0.5)");
        }else if (Boufe_avant > bouffe) {
            // ON est pas dans l'optimum, on essaye de bouger :
            // On bouge de i_mvt points de grille (à definir en fonction de resolution et vitesse des poissons)
            // RANDOM WalK : 
            //System.out.println("CE POISSON RANDOM WALK car le phyto diminue. Boufe_avant = " + Boufe_avant + " bouffe maintenant = " + bouffe + " ancienne dir_banc : " + dir_banc);
            dir_banc = nbre_au_hazard_entre_0_et_x(3);
            vit_banc_Bls = 0.5f;
        }
            else {
                      //  System.out.println("CE POISSON est BIEN; il continue doucement dans cette direction : " + dir_banc);
             vit_banc_Bls = 0.1f;
            }
        float km_par_dt_adultes = (float) (vit_banc_Bls*Body_length/100* dt_adultes * 3600);
        float nbgridcell_par_dt_adultes = (float) km_par_dt_adultes / (float) Population.mean_grid_size;

// ON TOURNE SI YA LA COTE DEVANT :
        boolean cote_devant = false;
        boolean cote_sur_chemin =false;
        // 1 - Ya-t-il une cote devant, dans la direction "dir_banc"?
int i;
float route;
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


        // TEST
/*
        if (Dataset_EVOL.isInWater(x,y+1)) {y +=1;
        //            System.out.println(" TEST 1 :  Dataset_EVOL.isInWater(x,y+1) = " + Dataset_EVOL.isInWater(x,y+1));
        }
        if (Dataset_EVOL.isInWater(x+1,y)) {x +=1;
        //              System.out.println(" TEST 2 :  Dataset_EVOL.isInWater(x+1,y) = " + Dataset_EVOL.isInWater(x+1,y));
        }

        // FIN TEST
        if ((!Dataset_EVOL.isInWater(x,y+1))&& (!Dataset_EVOL.isInWater(x+1,y)))
        {
        //               System.out.println("CE POISSON est a la côte = " + x + " ; y = " + y);
        //               System.out.println(" Dataset_EVOL.isInWater(x,y+1) = " + Dataset_EVOL.isInWater(x,y+1) + " ; Dataset_EVOL.isInWater(x+1,y) = " +  Dataset_EVOL.isInWater(x+1,y));
        }
         */
         

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


        if (!Dataset_EVOL.isInWater(x, y)){
            System.out.println("MORT A TERRE BUGBUGBUGBUGBUG DIR_BANC = " + dir_banc + "  ; nbgridcell_par_jour = " + nbgridcell_par_dt_adultes);
        }






