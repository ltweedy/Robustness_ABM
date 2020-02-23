/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:46
 * To change this template use File | Settings | File Templates.
 */

import java.util.ArrayList;
import java.util.Random;

public class Cell {
    public double[] position;
    public double[] force;
    public double[] oldforce;
    public double[] oldoldforce;
    public MigrationSimulation sim;
    public double meanOccupancy = 0.5;

    //public static double minSurvivalOccupancy = 0.25;
    public static double minMitogenicOccupancy = 0.35;
    public static double speed = 12;  // um/min
    public static double rspeed = Math.sqrt(speed);  // um/min
    public static double kBIRTH = 0.00;//0.05
    public static double kDEATH = 0.08;
    public static double minuteAverageOccupancy = 0.1*ChemicalEnvironment.sMax;

    long iTicks;

    public static double kD      = 0.05;
    public static double kI      = 0.02;
    public static int nReceptors = 10000000;
    public static Random RG = new Random();
    public boolean original = true;
    public DoubleRectangle box;

    public double minrad = Math.random()>0.00? 24.0 : 48;     // Below which pushing out
    public double neuralrad = 8.0;  // Above which pulling in
    public double maxrad = 10.0;  // Maximum, beyond which no forces


    public static boolean inducible = false;
    public boolean weakDegrader = false;

    public static double degR = 0.2; // Basal degradation
    public double deg = 0.0;


    public static double kRepel = 0.5;
    public static double kAttract = 0.1;

    public double ld  = 10.0;
    public static double CIbMax = 4;
    public double CIb = CIbMax;
    public double oF = 0;
    public double oB = 0;
    public double oU = 0;
    public double oL = 0;
    public boolean active = true;
    public boolean sticky = true;

    public boolean dgr = true;
    public ArrayList<EnvironmentPoint> points = new ArrayList<EnvironmentPoint>();




    //public double[] velocity;

    public double x(){ return position[0];}
    public double y(){ return position[1];}
    public double fx(){ return force[0];}
    public double fy(){ return force[1];}
    public double f2x(){ return oldforce[0];}
    public double f2y(){ return oldforce[1];}
    public double f3x(){ return oldoldforce[0];}
    public double f3y(){ return oldoldforce[1];}



    public Cell(double[] position, MigrationSimulation sim){

        this.position = position;
        this.force = new double[]{0.0,0.0};
        this.oldforce = new double[]{0.0,0.0};
        this.oldoldforce = new double[]{0.0,0.0};
        this.box = new DoubleRectangle(position[0]-maxrad,
                                       position[1]-maxrad,
                2.0*maxrad,
                2.0*maxrad   //width and height the same.
                                );

        this.sim = sim;
        this.original = true;
    }

    public void clear(){
        iTicks++;

        if(iTicks%200==0) {

            oldoldforce[0] = oldforce[0];
            oldoldforce[1] = oldforce[1];

            oldforce[0] = force[0];
            oldforce[1] = force[1];

        }
        force[0] = 0.0;
        force[1] = 0.0;

    }

    /*public boolean CheckSurvival(){
        if(meanOccupancy<minSurvivalOccupancy && Math.random()*AgentBasedSimulation.dt<kDEATH) return false;
        else return true;
    }
    */
    public boolean CheckMitosis(){

        double mitosisChance = (minuteAverageOccupancy-minMitogenicOccupancy)/(1-minMitogenicOccupancy);
        if(mitosisChance>0 && Math.random()<AgentBasedSimulation.dt*kBIRTH*mitosisChance){
            return true;
        }
        return false;
    }

    public void Refactory(){
        minuteAverageOccupancy = 0;
    }

    public void updatePosition(){

        if(!active) return;

        double npx = this.position[0]+force[0]* AgentBasedSimulation.dt;
        double npy = this.position[1]+force[1]* AgentBasedSimulation.dt;

        if(!sim.environment.GetCellPass(npx,npy)){

            force[0] *= -0.3;
            force[1] *= -0.3;

        }

        this.position[0]+=force[0]* AgentBasedSimulation.dt;
        this.position[1]+=force[1]* AgentBasedSimulation.dt;

        this.box.x = (position[0]-minrad);
        this.box.y = (position[1]-minrad);

        //if(this.position[1]<MigrationSimulation.padding||this.position[1]> sim.yMax-MigrationSimulation.padding) {
        //    this.position[1] = MyMaths.bounded(MigrationSimulation.padding, sim.yMax - MigrationSimulation.padding, this.position[1]);
        //}
        //if(this.position[1]<MigrationSimulation.padding) this.position[1] = 2.0*MigrationSimulation.padding-this.position[1];
        //else if(this.position[1]> sim.yMax-MigrationSimulation.padding) this.position[1] = 2.0* (sim.yMax-MigrationSimulation.padding) - this.position[1];


        //if(this.position[0]<MigrationSimulation.padding) this.position[0] = 2.0*MigrationSimulation.padding-this.position[0];
        //else if(this.position[0]> sim.xMax-MigrationSimulation.padding){
        //    this.position[0] = 2.0* (sim.xMax-MigrationSimulation.padding) - this.position[0];
            //System.out.println((sim.Ttotal/60));
        //}

    }

    public void addForce(double dx, double dy){
        force[0]+=dx;
        force[1]+=dy;
    }

    public void GetEnvironmentPointsInCell(ChemicalEnvironment env){
        points.clear();

        double dx = ChemicalEnvironment.grain;
        int iX = Math.max(0,(int) Math.floor((position[0]-minrad)/dx));
        int iY = Math.max(0,(int) Math.floor((position[1]-minrad)/dx));

        int iW = (int) Math.ceil(2.0*minrad/dx);

        for(int i = iX; i<=iX+iW; i++){
            for(int j = iY; j<=iY+iW; j++){
                if((((i*dx-position[0])*(i*dx-position[0])) + ((j*dx-position[1])*(j*dx-position[1])))<=minrad*minrad){
                    int i2 = Math.max(0,Math.min(env.profile.size()-1, i));
                    int j2 = Math.max(0,Math.min(env.profile.get(i2).size()-1, j));
                    EnvironmentPoint ep = env.profile.get(i2).get(j2);
                    points.add(ep);
                }
            }
        }

    }


    public void DegradeFromEnvironment(ChemicalEnvironment env){
        points.clear();
        GetEnvironmentPointsInCell(env);

        minuteAverageOccupancy = (1-AgentBasedSimulation.dt)*minuteAverageOccupancy + AgentBasedSimulation.dt*meanOccupancy;

        //double newS = ChemicalEnvironment.sMax*(0.05+0.95*minuteAverageOccupancy*minuteAverageOccupancy/(0.4*0.4+minuteAverageOccupancy*minuteAverageOccupancy));
        //double newS = ChemicalEnvironment.sMax*(0.05+0.95*minuteAverageOccupancy);
        //double newS = ChemicalEnvironment.sMax;
        double newS = 0;
        
        // How much of the cell's degrading power is in each point.
        double fr = newS/points.size();
        double km = ChemicalEnvironment.kM;
        double dt = AgentBasedSimulation.dt;


        double k1, k2, k3, k4;

        for(EnvironmentPoint e : points) {
            // If the point allows free, unrestricted diffusion, degrade.

            //RK4 integration.
            k1 = fr * e.c / (e.c + km);

            k2 = (e.c + 0.5 * dt * k1);
            k2 = fr * k2 / (k2 + km);

            k3 = (e.c + 0.5 * dt * k2);
            k3 = fr * k3 / (k3 + km);

            k4 = e.c + dt * k3;
            k4 = fr * k4 / (k4 + km);

            double degraded = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1 + 2 * k2 + 2 * k3 + k4);
            double d2 = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1 + 2 * k2 + 2 * k3);
            double d3 = /*minrad*minrad*minrad*0.5**/(dt / 6.0) * deg * (k1);


            if (!e.fixed && e.open) {
                e.c -= degraded;
                e.c_m1-=d2;
                e.c_m2-=d3;

                //e.c2 += degraded;

            }
        }
    }

    public void AddToEnvironment(ChemicalEnvironment env){

        for(EnvironmentPoint e : points) {
            // If the point allows free, unrestricted diffusion, degrade.
            e.c += 0.1*AgentBasedSimulation.dt;
            //e.c2 += 0.085*AgentBasedSimulation.dt;

        }
    }

    public double[] EstimateGradientDirection(ChemicalEnvironment env){

        GetEnvironmentPointsInCell(env);

        double estX = 0;
        double estY = 0;
        double cX = 0;
        double cY = 0;

        double xMin = 0;
        double xMax = 0;
        double yMin = 0;
        double yMax = 0;


        for(EnvironmentPoint ep : points){
            cX += ep.x;
            cY += ep.y;
        }

        cX/=points.size();
        cY/=points.size();

        double occupancy;
        meanOccupancy = 0;
        double nT;
        double nR;
        double dC;

        for(EnvironmentPoint ep : points) {
            // Calculate directional projection of occupancy across the cell
            occupancy           = (ep.c / (ep.c + kD/*(1+ep.c2/kI)*/));
            meanOccupancy      += occupancy/points.size();
            estX               += occupancy * (ep.x - cX) / points.size();
            estY               += occupancy * (ep.y - cY) / points.size();
        }

        if(inducible) {

            double x = meanOccupancy;
            //x = kD * x / (1-x);                                 // These run if
            //x = (x*x)/(kD+(x*x));                               // ultrasensitive.

            deg = x;                                            // If no delay
            //deg += degR * (x - deg) * MelaMigration.dt;       // If delay
        }
        else deg = 1;

        // Remove environmental attractant based on occupancy, area, depth, dt.

        //number of molecules in area
        /*    if(ep.open&&!ep.fixed) {

                nT = ChemicalEnvironment.environment_depth * ChemicalEnvironment.grain
                        * ChemicalEnvironment.grain * 1E-15 * ChemicalEnvironment.AVAGADRO * ep.c * 1E-6;
                // Number of occupied receptors in THIS part of the cell.
                nR = (nReceptors / points.size()) * occupancy;

                dC = Math.min(1, nR / nT) * ep.c * MelaMigration.dt;

                //System.out.println("nT -> "+nT+",    nR-> "+nR+",   dC -> "+dC);
                // Removed fraction of total concentration involved in gathering this information!
                ep.c -= dC;
            } */

        return new double[] {estX,estY};


        /*double cMax = 0;
        EnvironmentPoint epo = points.get(0);

        for(EnvironmentPoint ep : points){
            double cr = ep.c+Math.sqrt(ep.c)*RG.nextGaussian();

            if(cr>cMax){
                cMax = cr;
                epo = ep;
            }
        }

        if(cMax<0.01) return new double[]{0.,0.};
        else return new double[]{epo.x-this.x(),epo.y-this.y()};   */
    }

    public boolean checkExit(){
        for(EnvironmentPoint ep : points){
            if(ep.fixed) return true;
        }
        return false;
    }

    public void updateGrowth(){

        if(sim.dieoff)      this.cull(0.0,0.0,4.0);
        if(sim.proliferate) this.split(0.0,0.0, AgentBasedSimulation.dimensions[0]);

    }


    public synchronized void cull(double a, double b, double c){

        double t1 = a*x()*x();
        double t2 = b*x();
        double t3 = c;

        if((kDEATH* AgentBasedSimulation.dt>Math.random())&&((t1+t2+t3)<Math.random()* AgentBasedSimulation.dimensions[0])){
            sim.cells.remove(this);
        }
    }

    public synchronized void split(double a, double b, double c){

        if((kBIRTH* AgentBasedSimulation.dt>Math.random())/*&&((t1+t2+t3)>Math.random()*MelaMigration.dimensions[0])*/){
            Cell clone = new Cell(new double[] {position[0] + 2.0*minrad*RG.nextGaussian(), position[1] + 2.0*minrad*RG.nextGaussian()},this.sim);
            clone.position[1] = MyMaths.bounded(0.0, AgentBasedSimulation.dimensions[1]-0.0, clone.position[1]);
            if(clone.position[0]<0.1) clone.position[0] = 0.2-clone.position[0];
            else if(clone.position[0]> AgentBasedSimulation.dimensions[0]-0.1) clone.position[0] = 2.0* (AgentBasedSimulation.dimensions[0]-0.1) - clone.position[0];


            sim.newCells.add(clone);
            clone.original = false;
        }
    }
}
