import java.io.File;
import java.io.IOException;
import java.util.Random;


/**
 * Created with IntelliJ IDEA.
 * User: luke
 * Date: 18/02/2014
 * Time: 08:41
 * To change this template use File | Settings | File Templates.
 */
public class AgentBasedSimulation
{

    //public static double  dimensions[] = {1000, 400, 1000, 400}; // um
    //public static double  dimensions[] = {3500, 3500, 3500, 3500}; // um
    public static double  dimensions[] = {6000, 400, 6000, 400}; // um
    //public static double  dimensions[] = {2000, 1800, 2000, 1800}; // um
    public static boolean visualise = true;
    public static boolean record = true;
    public static boolean exdeg  = true;

    public static int     yPos = 40;
    public static boolean abs = true;
    public static boolean pinned = true;
    public static boolean ctc = false;
    public static boolean competivite = false;

    static int pop = 800;                    // Initial population
    static double T  = 7.51*1.0*60;                // Days, Hours, Minutes
    public static double tTotal = 0;
    public static double dt = 0.05; //0.0015        // 0.05~30.0s
    public static double alpha = 0.6;
    public static double outInt = 0.5;        // Every 15 sec.
    public static double rdt = Math.sqrt(dt);   // Sqrt of dt for brownian motion
    public static String directory = System.getProperty("user.home")+"/Science/RobustnessConcentration/LongNoDegExample/10uM2/";
    public static boolean finished = false;

    public Thread visualiser;
    //public ExecutorService threadpool;
    //ArrayList<Callable<Object>> tasks;

    public MigrationSimulation RWS;

    private static void RandomiseMazeProperties(){

        Random r = new Random();

        double exp = 1+r.nextDouble()*6.5;
        Cell.speed = Math.pow(exp, 2);    // 0.04-100 spaced geometrically.
        ChemicalEnvironment.L = 100*(1+r.nextInt(20)); // 100-2000
        ChemicalEnvironment.W = 10  + r.nextInt(66);    // 10-75
        ChemicalEnvironment.DiffC = Math.pow(50 + r.nextInt(120), 2);
        ChemicalEnvironment.baseConcentration = Math.pow(10, -1 + r.nextDouble() * 3.5);  // 0.1 to 1000 exponentially
        AgentBasedSimulation.dt = 2.5E-0/ ChemicalEnvironment.DiffC; //
        AgentBasedSimulation.dt = Math.min(AgentBasedSimulation.dt, 5E-2);
        System.out.println("L-> " + ChemicalEnvironment.L + ",  W-> " + ChemicalEnvironment.W +
                ",  D-> " + ChemicalEnvironment.DiffC + ",  C0-> " + ChemicalEnvironment.baseConcentration +
                ",  S-> " + Cell.speed + ", dt-> " + AgentBasedSimulation.dt);
    }

    private AgentBasedSimulation(boolean skipguis){

        if(record){
             boolean bDir = new File(directory).mkdirs();
        }
        /*if(visualise && !skipguis){
           makeSimGUIs();
        }
        else{*/
            RWS = new MigrationSimulation(false,false,ctc,abs,alpha, Cell.speed,dt,ChemicalEnvironment.grain, ChemicalEnvironment.DiffC, Cell.kD, ChemicalEnvironment.kM, ChemicalEnvironment.sMax);
        //}
        Thread th1 = new Thread() {
            public void run() {
                long startTime = 0;
                long endTime   = 0;

                int j = 0;
                for (int i = 0; i <= T / dt; i++) {
                    //if(i==0)    startTime = System.currentTimeMillis();
                    /*else if(i%5000 == 0) {
                        endTime = System.currentTimeMillis() - startTime;
                        System.out.println("TIME 1:: "+endTime);
                        startTime = System.currentTimeMillis();
                    }*/


                    if(RWS.complete) break;
                    if (RWS.paused) {
                        i--;
                    }
                    else{
                        RWS.Ttotal += dt;

                        j++;
                        RWS.Iterate();
                        if (i % 1000 == 1 && visualise) RWS.controlPanel.t_time.setText(Double.toString(dt * i / 60.0).substring(0,3));
                        if (j * dt >= outInt) {
                            RWS.paused = true;
                            j = 0;
                            RWS.SnapImage("con");
                            RWS.WriteCellData();
                            RWS.WriteEnvironmentData();
                            RWS.paused = false;
                        }

                    }
                }

                if(record) {
                    RWS.WriteCellData();
                    RWS.WriteEnvironmentData();
                    RWS.WriteShortRecord();
                    try {
                        RWS.bw.flush();
                        RWS.bw.close();
                        RWS.fw.close();
                        RWS.bw4.flush();
                        RWS.bw4.close();
                        RWS.fw4.close();
                        RWS.bw3.flush();
                        RWS.bw3.close();
                        RWS.fw3.close();
                    } catch (IOException e) {
                    }
                }
                finished = true;
                System.exit(1);
            }
            //if(visualise) while(true)    try{sleep(100);} catch(Exception e){}
        };
        th1.run();
    }

    public static void main(String[] args){

        if(args.length==0) new AgentBasedSimulation(false);
        if(args.length==1){
            if(args[0].toString().equals("R")) {
                RandomiseMazeProperties();
                //visualise = false;
                new AgentBasedSimulation(true);
            }
        }
        else if(args.length%2==0){
            parseArgs(args);
            visualise = false;
            new AgentBasedSimulation(false);
        }
        else System.exit(-1);
    }

    private void makeSimGUIs(){

        SimGUIPanel s1 = new SimGUIPanel();

        s1.create();

        RWS = s1.SetupSimulation("/1.txt");

    }

    private static void parseArgs(String[] args){

        for(int i=0; i<args.length; i+=2){
            try{
                if(args[i].equals("p"))             pop = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("D"))        ChemicalEnvironment.DiffC = Double.parseDouble(args[i+1]);
                else if(args[i].equals("kD"))       Cell.kD    = Double.parseDouble(args[i+1]);
                else if(args[i].equals("kM"))       ChemicalEnvironment.kM    = Double.parseDouble(args[i+1]);
                else if(args[i].equals("sMax"))     ChemicalEnvironment.sMax  = Double.parseDouble(args[i+1]);
                else if(args[i].equals("a"))        abs = Boolean.parseBoolean(args[i+1]);
                else if(args[i].equals("c"))        ctc = Boolean.parseBoolean(args[i+1]);
                else if(args[i].equals("in"))       MigrationSimulation.sMazePicture = args[i+1];
                else if(args[i].equals("out"))      directory = args[i+1];
                else if(args[i].equals("alpha"))    alpha = Double.parseDouble(args[i+1]);
                else if(args[i].equals("k"))        ChemicalEnvironment.k = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("L"))        ChemicalEnvironment.L = Double.parseDouble(args[i + 1]);
                else if(args[i].equals("skew"))     ChemicalEnvironment.skew = Integer.parseInt(args[i + 1]);
                else if(args[i].equals("asym"))     ChemicalEnvironment.asym = Double.parseDouble(args[i + 1]);
            }
            catch(NumberFormatException e){}
        }

    }
}



