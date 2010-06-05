package edu.mit.csail.cgs.deepseq.features;

import java.util.HashMap;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.StrandedBase;

/**
 * BindingComponents are used in mixture models to represent potential binding events.
 * @author shaun
 *
 */
public class BindingComponent implements Comparable<BindingComponent>{

	private BindingModel model;
	private Point position;
	private double mixingProb;
	private double[] conditionBeta;
	private int numConditions;
	private HashMap<StrandedBase,Double>[] responsibilities;
	private double[][] readProfile_plus;
	private double[][] readProfile_minus;
	private Point EM_position;		// EM result
	private int old_index;
	private double alpha;

	public int getOld_index() {
		return old_index;
	}
	public void setOld_index(int old_index) {	// store the original index as initialized
		this.old_index = old_index;
	}
	public Point getEMPosition() {
		if (EM_position==null)
			return position;
		else
			return EM_position;
	}
	public void setEMPosition(Point EM_position) {
		this.EM_position = EM_position;
	}
	public BindingComponent(BindingModel m, Point pos){this(m,pos,1);}
	public BindingComponent(BindingModel m, Point pos, int numConds){
		model=m;
		position=pos;
		numConditions = numConds;
		conditionBeta=new double[numConditions];
		// uniform initialization
		double uc=1/(double)numConditions;
		for(int c=0; c<numConditions; c++){
			conditionBeta[c]=uc;
		}
		mixingProb=1;
		responsibilities = new HashMap[numConditions];
		for(int c = 0; c < numConditions; c++) { responsibilities[c] = new HashMap<StrandedBase,Double>(); /*normalized*/ }			
	}
	
	//Accessors
	public double scoreHit(ReadHit h){
		int dist = h.getStrand()=='+' ? h.getFivePrime()-position.getLocation():position.getLocation()-h.getFivePrime();
		return model.probability(dist);
	}
	public double scoreBase(StrandedBase b){
		int dist = b.getStrand()=='+' ? b.getCoordinate()-position.getLocation():position.getLocation()-b.getCoordinate();
		return model.probability(dist);
	}
	public double score(int dist){
		return model.probability(dist);
	}
	public Point getLocation(){return position;}
	public int getNumConditions(){return numConditions;}
	public double getMixProb(){return mixingProb;}

	public double getSumResponsibility(int cond){
		double sum = 0;
		for (StrandedBase base:responsibilities[cond].keySet())
			sum += base.getCount()*responsibilities[cond].get(base);
		return(sum);
	}
	
	public double getTotalSumResponsibility(){
		double sum = 0;
		for(int c = 0; c < numConditions; c++) { sum += getSumResponsibility(c); }
		return sum;
	}
	

	public double getConditionBeta(int c){
		if(c>=numConditions){return -1;}
		else{
			return(conditionBeta[c]);
		}
	}
	public double[][] getReadProfile_plus(){
		return readProfile_plus;
	}
	public double[][] getReadProfile_minus(){
		return readProfile_minus;
	}
	public double[] getReadProfile_plus(int condition){
		return readProfile_plus[condition];
	}
	public double[] getReadProfile_minus(int condition){
		return readProfile_minus[condition];
	}	
	public double getReadProfile(int cond, int index, char strand){
		double result=0;
		if (strand=='+')
			result=readProfile_plus[cond][index];
		else if (strand=='-')
			result=readProfile_minus[cond][index];
		return result;
	}
	//Mutators
	public void setMixProb(double p){mixingProb=p;}
	public void setResponsibility(int cond, StrandedBase b, double x){responsibilities[cond].put(b, x);}
	public void setConditionBeta(int cond, double beta){
		if(cond<numConditions){
			conditionBeta[cond]=beta;
		}
	}
	public void setReadProfile(int cond, double[] profile, char strand){
		if (readProfile_plus==null){
			readProfile_plus = new double[numConditions][profile.length];
			readProfile_minus = new double[numConditions][profile.length];
		}
		if (strand=='+')
			readProfile_plus[cond] = profile;
		else if (strand=='-')
			readProfile_minus[cond] = profile;
	}

	//Comparable default method
	public int compareTo(BindingComponent m) {
//		return compareByAvgQValue(f);
		return getLocation().compareTo(m.getLocation());
	}
	public void setAlpha(double alpha) {
		this.alpha = alpha;		
	}
	public double getAlpha(){
		return alpha;
	}
	public String toString(){
		return getLocation().getLocationString()+" EM:"+getEMPosition().getLocationString();
	}
}