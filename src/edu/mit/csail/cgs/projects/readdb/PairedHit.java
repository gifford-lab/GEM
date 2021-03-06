package edu.mit.csail.cgs.projects.readdb;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class PairedHit implements Comparable<PairedHit> {

    public int leftChrom, rightChrom;
    public int leftPos, rightPos;
    public float weight;
    public boolean leftStrand, rightStrand;
    public short leftLength, rightLength;

    public PairedHit(int leftchrom, int leftpos, boolean leftstrand, short leftlen, 
                     int rightchrom, int rightpos, boolean rightstrand, short rightlen,
                     float weight) {
        leftChrom = leftchrom;
        rightChrom = rightchrom;
        leftPos = leftpos;
        rightPos = rightpos;
        this.weight = weight;
        leftStrand = leftstrand;
        rightStrand = rightstrand;
        leftLength = leftlen;
        rightLength = rightlen;            
    }

    public boolean equals(Object o) {
        if (o instanceof PairedHit) {
            PairedHit other = (PairedHit)o;
            return (leftChrom == other.leftChrom &&
                    rightChrom == other.rightChrom &&
                    leftPos == other.leftPos &&
                    rightPos == other.rightPos &&
                    leftStrand == other.leftStrand &&
                    rightStrand == other.rightStrand &&
                    Math.abs(weight - other.weight) < .001);
        } else {
            return false;
        }
    }   
    
    @Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + leftChrom;
		result = prime * result + leftPos;
		result = prime * result + (leftStrand ? 1231 : 1237);
		result = prime * result + rightChrom;
		result = prime * result + rightPos;
		result = prime * result + (rightStrand ? 1231 : 1237);
		return result;
	}
    //Sort using left only
    public int compareTo(PairedHit o) {
        if (leftChrom == o.leftChrom) {
            return leftPos - o.leftPos;
        } else {
            return leftChrom - o.leftChrom;
        }
    }
    
    public StrandedRegion leftRegion(Genome genome) {
    	return new StrandedRegion(genome, genome.getChromName(leftChrom), leftPos, leftPos+(leftStrand ? leftLength : -leftLength), leftStrand ? '+' : '-');
    }
    
    public StrandedRegion leftPointRegion(Genome genome) {
    	return new StrandedRegion(genome, genome.getChromName(leftChrom), leftPos, leftPos+1, leftStrand ? '+' : '-');
    }
    
    public StrandedPoint leftStrandedPoint(Genome genome) {
    	return new StrandedPoint(genome, genome.getChromName(leftChrom), leftPos, leftStrand ? '+' : '-');
    }
    
    public StrandedRegion rightRegion(Genome genome) {
    	return new StrandedRegion(genome, genome.getChromName(rightChrom), rightPos, rightPos+(rightStrand ? rightLength : -rightLength), rightStrand ? '+' : '-');
    }
    
    public StrandedRegion rightPointRegion(Genome genome) {
    	return new StrandedRegion(genome, genome.getChromName(rightChrom), rightPos, rightPos+1, rightStrand ? '+' : '-');
    }
    
    public StrandedPoint rightStrandedPoint(Genome genome) {
    	return new StrandedPoint(genome, genome.getChromName(rightChrom), rightPos, rightStrand ? '+' : '-');
    }
    
    public boolean leftContainedIn(Region r) {
    	return r.contains(leftPointRegion(r.getGenome()));
    }
    
    public boolean rightContainedIn(Region r) {
    	return r.contains(rightPointRegion(r.getGenome()));
    }
    
    public int lesserPos() {
    	if (leftChrom == rightChrom) {
    		return leftPos < rightPos ? leftPos : rightPos;
    	} else {
    		return leftChrom < rightChrom ? leftPos : rightPos;
    	}
    }
    
    public int greaterPos() {
    	if (leftChrom == rightChrom) {
    		return leftPos > rightPos ? leftPos : rightPos;
    	} else {
    		return leftChrom > rightChrom ? leftPos : rightPos;
    	}
    }
    
    public boolean lesserStrand() {
    	if (leftChrom == rightChrom) {
    		return leftPos < rightPos ? leftStrand : rightStrand;
    	} else {
    		return leftChrom < rightChrom ? leftStrand : rightStrand;
    	}
    }
    
    public boolean greaterStrand() {
    	if (leftChrom == rightChrom) {
    		return leftPos > rightPos ? leftStrand : rightStrand;
    	} else {
    		return leftChrom > rightChrom ? leftStrand : rightStrand;
    	}
    }
    
    public PairedHit flippedCopy() {
    	return new PairedHit(rightChrom, rightPos, rightStrand, rightLength, leftChrom, leftPos, leftStrand, leftLength, weight);
    }

	public void flipSides() {
        int x = leftChrom;
        leftChrom = rightChrom;
        rightChrom = x;

        x = leftPos;
        leftPos = rightPos;
        rightPos = x;

        boolean b = leftStrand;
        leftStrand = rightStrand;
        rightStrand = b;

        short s = leftLength;
        leftLength = rightLength;
        rightLength = s;
    }
    public String toString() {
        return String.format("%d:%d,%d:%c and %d:%d,%d:%c weight %.2f",
                             leftChrom, leftPos, leftLength, leftStrand ? '+' : '-',
                             rightChrom, rightPos, rightLength, rightStrand ? '+' : '-',
                             weight);
    }

}