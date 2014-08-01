/*
 * Copyright (C) 2012 tbrochier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package evolution;

//import evolution.Poissons;
/**
 *
 * @author tbrochier,lpecquerie
 */
public class DebParticleLayer {

    //  double dt; // pas de temps
    /**
     * Particle length [millimeter]
     */
    private float length;
    private double weigth;
    public int stade, Nb_oeuf; // 0 = EGG; 1 = YOLK_SAC_LARVA; 2 = FEEDING_LARVA; 3 = JUVENILE;
    /*
     * DEB VARIABLES AND PARAMETERS
     */
    //DEB VARIABLES
    double E_init;
    public double E; // Reserve
    public double V; // Structure
    public double E_H; // Maturity
    public double E_R; // Reproduction
    // [TIM]
    public double Gonade_buffer; // Dynamique de la gonade
    public double f;
    // DEB PARAMETERS
    // temperature correction
    static double T_ref;    // K, Reference temperature ;
    static double T_A;      // K, Arrhenius temperature ;
    static double T_L;     // K, lower boundary tolerance range
    static double T_H;     //  K, upper boundary tolerance range
    static double T_AL;    //  K, Arrhenius temp for lower boundary
    static double T_AH;   //   K, Arrhenius temp for upper boundary
    // feeding & assimilation
    double F_m;    //  l/d.cm^2, {F_m} max spec searching rate
    double kap_X;         //  -, digestion efficiency of food to reserve
    double p_Am;     //  J/cm^2/d, maximum surface-specific assimilation rate
    // mobilisation, maintenance, growth & reproduction
    double v;        //  cm/d, energy conductance
    double Kappa;      //  -, allocation fraction to soma = growth + somatic maintenance
    double kap_R;    //  -, reproduction efficiency
    double p_M;      //  J/d.cm^3, [p_M], vol-specific somatic maintenance
    double p_T;      //  J/d.cm^2, {p_T}, surface-specific som maintenance
    double k_J;      //  1/d, maturity maint rate coefficient
    double E_G;      //  J/cm^3, [E_G], spec cost for structure
    // life stages: E_H is the cumulated energy from reserve invested in maturation
    double E_Hh;	// J, Maturity threshold at hatching
    double E_Hb;        // J, E_H^b, maturity at birth
    double E_Hj;	// J, Maturity threshold at metamorphosis
    double E_Hp;       // J, E_H^p, maturity at puberty
// [TIM]
    double E_Hg;      //J, E_H^g, maturité des gonades
    double E_Hbatch;      //J, E_H^g, maturité d'un batch de ponte
    // param to compute observable quantities
    double del_M;  //  -, shape coefficient to convert vol-length to physical length
    double d_V;  //  g/cm^3, specific density of structure (dry weight)
    double mu_V;  //  J/mol, specific chemical potential of structure
    double mu_E;  //  J/mol, specific chemical potential of reserve
    double w_V;  //  g/mol, molecular weight of structure
    double w_E;  //  g/mol, molecular weight of reserve
    double c_w;  //  -, water content (c_w * W_w = total water weight)
    // compounds parameters
    double X_K;  // same unit as food, half-saturation coefficient
    double p_Xm; // J.cm-2.d-1, max surf area specific ingestion rate, here we assume constant assimilation efficiency
    // assimilation parameter = primary parameter, ingestion rate is calculated "backward" from assimilation
    double W_V, W_E, W_ER; // poids en gramme de chaque compartiment, pour le calcult du poids total
    // Tim parameters
    boolean adult_ready_to_spawn;
    int Nb_eggs_min, relative_fecondity; // le nombre minimal d'oeuf par evenement de ponte (proportionnel a la taille du poisson??)

    public void init() {
        loadParameters();
        // CONDITIONS INITIALES
        E_init = 2.58; // J, Initial Reserve = egg
        E = E_init; // J, Initial Reserve = egg
        V = 0.000001; // cm^3, Initial volume --> close to 0, try different initial values
        // taille d'un oeuf = 1 mm3
        E_H = 0; //
        E_R = 0;

        //[TIM]
        Gonade_buffer = 0;

    }

    private void loadParameters() {
        //Primary parameters
        // temperature correction
        T_ref = 20 + 273;   //  K, Reference temperature ; 
        T_A = 5000;//9800;       //  K, Arrhenius temperature ;
        T_AL = 50000;      //  K, Arrhenius temp for lower boundary
        T_AH = 190000;     //  K, Arrhenius temp for upper boundary

// GAMME MIN-MAX DE TEMP ICI CELLE DES ADULTES INDICATIF,
        //MAIS CHANGE POUR CHAQUE STADE (dans set_stage):
        // --> POUR CHANGER CES VALEURS LE FAIRE DANS set_stage
//        T_L = 10 + 273;//0 + 273;    //  K, lower boundary tolerance range
  //      T_H = 29 + 273;//27 + 273;   //  K, upper boundary tolerance range


        // feeding & assimilation
        F_m = 6.51;    // l/d.cm^2, {F_m} max spec searching rate
        kap_X = 0.8;   // -, digestion efficiency of food to reserve
        p_Am = 1.677 * 92.51 / 0.3436 * 2.4019;     // J/cm^2/d, maximum surface-specific assimilation rate

        // mobilisation, maintenance, growth & reproduction
        v = 0.1379 * 2.4019;  // cm/d, energy conductance
        Kappa = 0.3436;       // -, allocation fraction to soma = growth + somatic maintenance
        kap_R = 0.95;         // -, reproduction efficiency
        p_M = 92.51;          // J/d.cm^3, [p_M], vol-specific somatic maintenance
        p_T = 0;              // J/d.cm^2, {p_T}, surface-specific som maintenance
        k_J = 0.002;          // 1/d, maturity maint rate coefficient
        E_G = 4767;          // J/cm^3, [E_G], spec cost for structure

        // life stages: E_H is the cumulated energy from reserve invested in maturation
        E_Hh = 1;        // J, E_H^h, maturity at hatching
        E_Hb = 1.372e0;  // J, E_H^b, maturity at birth
        E_Hj = 20;       // J, E_H^j, maturity at metamorphosis
        E_Hp = 1.928e5;  // J, E_H^p, maturity at puberty


        // essais E_Hp = 1.928e6 simu 39 --> pas de maturation avant 2 ans (deux ans de run)
        // essais E_Hp = 1.928e5 simu 39 --> pas de maturation avant 6 mois


        // param to compute observable quantities
        del_M = 0.1391;    //  -, shape coefficient to convert vol-length to physical length // maitenant variable selon le stade (voir setStage())
        d_V = 0.2; 	   // g/cm^3, specific density of structure (dry weight)
        mu_V = 500000;    // J/mol, specific chemical potential of structure
        mu_E = 550000;    // J/mol, specific chemical potential of reserve
        w_V = 23.9;      // g/mol, molecular weight of structure
        w_E = 23.9;     // g/mol, molecular weight of reserve
        c_w = 1 - d_V;   // - , water content (c_w * W_w = total water weight)
        relative_fecondity = 400; // nb oeuf par grammes de femelle These Freon

        // compound parameters
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        X_K = 3; // 0.2 = bon pour le Nano_phyto en surface //(p_Am / (kap_X * F_m))/100;// c'etait 50  // same unit as food, half-saturation coefficient
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p_Xm = p_Am / kap_X; // J.cm-2.d-1, max surf area specific ingestion rate, here we assume constant assimilation efficiency

        // time step
        //dt = 1; // 1 JOUR           Simulation.dt_advec; // AVERIFIER --- - - -

    }

    public void computeLength() { /// LAURE: shape a modifier
        length = (float) (Math.pow(V, 1 / 3.0) / del_M);
    }

    public float getLength() {
        return length;
    }

    public void setLength(float length) {
        this.length = length;
    }

    public void computeWeight() {
        W_V = d_V * V;
        W_E = (w_E / mu_E) * E;
        W_ER = (w_E / mu_E) * E_R;
        weigth = W_V + W_E + W_ER;
    }

    public double getWeigth() {
        return weigth;
    }

    public void setWeigth(double weigth) { // A VOIR A QUOI ça SERT??
        this.weigth = weigth;
    }

    public double getE() {
        return E;
    }

    public void setE(double E) {
        this.E = E;
    }

    public double getV() {
        return V;
    }

    public void setV(double V) {
        this.V = V;
    }

    public double getE_R() {
        return E_R;
    }

    public void setE_R(double E_R) {

        this.E_R = E_R;
    }

    public double getE_H() {
        return E_H;
    }

    public void setE_H(double E_H) {
        this.E_H = E_H;
    }

    public void setStage_properties() {
        // Dans cette routine on va mettre toutes les propriété qui changent
        // selon le stade (shape coefficient, ..?)
        // on peut dire qu'on appelle cette fonction tout les jour ou tout
        // les pas de temps (c juste un test --> ca coute rien)
        // je l'appelle donc ligne 226 avant computelength();
        // La variable publique stade pourra être acceder par poisson pour
        // ajuster son comportement mais il faudrais avoir 2 types d'adultes et juvenile (si les juvenile peuvent déja se reproduire?)
        // genre "adult_ready_to_spawn" et "adult_not_ready_to_spawn"
        // c'est possible? :-)

        /** Egg */
        if (E_H < E_Hh) {
            //    return Stage.EGG;
            del_M = 0.1391; // shape coefficient for EGG
            stade = 0;

            // Min et Max des corr de flux de temperature :
            T_L = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range

        } /** Yolk-sac-larva */
        else if ((E_Hh <= E_H) && (E_H < E_Hb)) {
            //    return Stage.YOLK_SAC_LARVA;
            del_M = 0.1391; // shape coefficient for YOLK_SAC_LARVA
            stade = 1;
            // Min et Max des corr de flux de temperature :
            T_L = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
        } /** Feeding larva */
        else if ((E_Hb <= E_H) && (E_H < E_Hj)) {
//            return Stage.FEEDING_LARVA;
            del_M = 0.1391;// shape coefficient for FEEDING_LARVA
            stade = 2;

            // Min et Max des corr de flux de temperature :
            T_L = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
        } /** Juvenile */
        else if ((E_Hj <= E_H) && (E_H < E_Hp)) {
            // return Stage.JUVENILE;
            del_M = 0.1391;// shape coefficient for JUVENILE
            stade = 3;
            // Min et Max des corr de flux de temperature :
            T_L = 15 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 27 + 273;//27 + 273;   //  K, upper boundary tolerance range
        } else {
            /** adult */
            //  return Stage.ADULT;
            del_M = 0.1391;// shape coefficient for ADULT
            stade = 4;
            // Min et Max des corr de flux de temperature :
            T_L = 10 + 273;//0 + 273;    //  K, lower boundary tolerance range
            T_H = 29 + 273;//27 + 273;   //  K, upper boundary tolerance range
        }

    }
// [TIM] : je rajoute les methodes de DebGrowthAction ici :

    public void execute(double temp, double food, double dt) {
        double[] res_deb = growth(dt, temp, food);
        E = res_deb[0];
        V = res_deb[1];
        E_H = res_deb[2];
        E_R = res_deb[3];

        setStage_properties();
        computeLength();
        computeWeight();


        if (stade == 4) {
            computeReproductionCycle();
        }

        if (res_deb[4] == 0) {
            //particle.kill(ParticleMortality.STARVATION);
        }
    }

    void computeReproductionCycle() {
        if (stade == 4) {

            Nb_eggs_min = (int) Math.round(weigth * relative_fecondity);

            // [TIM]
            // Quelle taille/energie de la gonade ?
            float alpha = 1f;
            E_Hg = alpha * E; //taille de la gonade proportionelle a la reserve (??)
            // E_Hbatch : énergie pour un batch de ponte de 400 oeufs/grammes :
            E_Hbatch = Nb_eggs_min * E_init / kap_R;

            if (E_Hg > 3 * E_Hbatch) { // si la gonade est potentielement assez grosse pour faire au moins 3 batch de ponte

             //   System.out.println("Gonade size = " + E_Hg + " ; Gonade_buffer = " + Gonade_buffer + " ; E_Hbatch = " + E_Hbatch + " ; E_R = " + E_R);

                if (Gonade_buffer < E_Hbatch) { // Si les gonades ne sont pas matures ou videe après une serie de ponte..
             //       System.out.println(" Remplissage gonade.. ");

                    if (E_R > E_Hg) { // si l'énergie le permet..
             //           System.out.println(" Gonade prête ");

                        Gonade_buffer = E_R; // Le Buffer de Gonade est rempli
                        E_R = E_R - E_Hg; // Et E_R est remis a 0
                    }
                }
            }

// sinon, si le buffer de Gonade EST remplis..
            // Les gonades sont prête pour un cycle de batch spawning
            // Prêt à pondre si il y a suffisament d'énergie pour un bacth

            adult_ready_to_spawn = ((E_R > E_Hbatch) && (Gonade_buffer >= E_Hbatch));

    //        System.out.println(" adult_ready_to_spawn = " + adult_ready_to_spawn);

        }
    }

    private double[] growth(double dt, double temperature, double food) {

        // if (V < Vj){ // no feeding < size-at-mouth opening
        if (E_H < E_Hb) { // no feeding < size-at-mouth opening
            f = 0;
        } else {
            f = food / (food + X_K);// scaled functional response
        }


        double cT = getcT(temperature);

        double p_XmT = p_Xm * cT;
        double p_AmT = p_XmT * kap_X;
        double v_T = v * cT;
        double p_MT = p_M * cT;
        double p_TT = p_T * cT;
        double k_JT = k_J * cT;


        //ENERGY FLUXES (J d-1)
        double flow_p_A = f * p_AmT * Math.pow(V, 2 / 3.0); // Assimilation
        double flow_p_M = p_MT * V; // Volume-linked somatic maintenance
        double flow_p_T = p_TT * Math.pow(V, 2 / 3.0); // Surface-linked somatic maintenance
        double flow_p_S = flow_p_M + flow_p_T; // Somatic maintenance
        double flow_p_C = (E / V) * (E_G * v_T * Math.pow(V, 2 / 3.0) + flow_p_S) / (Kappa * E / V + E_G); // eq. 2.12 p.37 Kooijman 2010; pC = [pC]*V
        double flow_p_G = Math.max((double) 0, Kappa * flow_p_C - flow_p_S);   // Growth
        double flow_p_J = k_JT * E_H; // Maturity maintenance
        double flow_p_R = ((1 - Kappa) * flow_p_C) - flow_p_J; // Maturation or reproduction

        //// STATE VARIABLES - Differential equations ////////////////////////////////////////
        double dEdt = flow_p_A - flow_p_C;   // Reserve, J d-1 ; dE/dt = pA - pC;
        double dVdt = flow_p_G / E_G;           // Structure, cm^3 d-1 ; dV/dt = (kap * p_C - p_S)/EG;
        double dE_Hdt, dE_Rdt;
        if (E_H < E_Hp) {
            dE_Hdt = flow_p_R; // maturation, J d-1
            dE_Rdt = 0; // Repro buffer J d-1
        } else {
            dE_Hdt = 0; // maturation, J d-1
            dE_Rdt = flow_p_R; // Repro buffer J d-1
        }

        //Integration
        E = E + dEdt * dt;
        V = V + dVdt * dt;
        E_H = E_H + dE_Hdt * dt;
        E_R = E_R + dE_Rdt * dt;

        // compute weight
        //double dV = 1;
        //W_dw[j] = V * dV  + (E/mu_E);
        //Compute DRY weight (g, dw) * 4.1 = Wet weight

        // starvation test
        int starvation;
        if (Kappa * flow_p_C < flow_p_M && (1 - Kappa) * flow_p_C < flow_p_J) {
            starvation = 0;//no starvation
        } else {
            starvation = 1;
        }

        double[] res = {E, V, E_H, E_R, (double) starvation};
        return res;
    }

    public void spawn() { // EN COURS
        Nb_oeuf = (int) Math.round(kap_R * E_R / E_init);
        Gonade_buffer = Gonade_buffer - E_Hbatch; // on entame la gonade
        E_R = E_R -E_Hbatch;
    }

    public static double getcT(double temperature) {
        // Correction of physiology parameters for temperature :
        double tempK = 273 + temperature;
        //double cT = Math.exp(T_A/T_ref-T_A/(tempK));
        double cT = Math.exp(T_A / T_ref - T_A / (tempK))
                * (1 + Math.exp(T_AL / T_ref - T_AL / T_L) + Math.exp(T_AH / T_H - T_AH / T_ref))
                / (1 + Math.exp(T_AL / tempK - T_AL / T_L) + Math.exp(T_AH / T_H - T_AH / tempK));
        return cT;
    }
}
