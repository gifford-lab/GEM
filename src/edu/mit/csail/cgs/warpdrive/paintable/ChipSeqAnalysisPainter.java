package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipSeqAnalysisModel;

public class ChipSeqAnalysisPainter extends RegionPaintable {


    private ChipSeqAnalysis analysis;
    private ChipSeqAnalysisModel model;
    private ChipSeqAnalysisProperties props;
    private DynamicAttribute attrib;
    private NonOverlappingLayout<ChipSeqAnalysisResult> layout;
    private ColorSet cs;

    public ChipSeqAnalysisPainter(ChipSeqAnalysis a, ChipSeqAnalysisModel m) {
        super();
        analysis = a;
        model = m;
        model.addEventListener(this);
        props = new ChipSeqAnalysisProperties();
        attrib = DynamicAttribute.getGlobalAttributes();
        layout = new NonOverlappingLayout<ChipSeqAnalysisResult>();
        cs = new ColorSet();
    }

    public ChipSeqAnalysisProperties getProperties() {return props;}
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        //        System.err.println("CCP.canpaint " + canPaint());
        if (!canPaint()) {
            return;
        }
        Collection<ChipSeqAnalysisResult> results = model.getResults();
        layout.setRegions(results);

        int numTracks = layout.getNumTracks();
        int w = lrx - ulx, h = lry - uly;
        int trackHeight = numTracks > 0 ? Math.max(1, h / numTracks) : 1;
        int spacing = Math.max(2, trackHeight/10);
        Region region = model.getRegion();
        int start = region.getStart(), end = region.getEnd();

        // set the min-width.
        double basesPerPixel = (double)(end-start) / (double)w;
        
        cs.reset();
        for (ChipSeqAnalysisResult r : results) {
            int x1 = getXPos(r.getStart(),
                             region.getStart(),
                             region.getEnd(),
                             ulx, lrx);
            int x2 = getXPos(r.getEnd(),
                             region.getStart(),
                             region.getEnd(),
                             ulx, lrx);
            int track = layout.getTrack(r);
            int y1 = uly + trackHeight * track;            
            g.setColor(cs.getColor());
            g.fillRect(x1,y1,x2-x1,trackHeight-spacing);
        }




    }

}