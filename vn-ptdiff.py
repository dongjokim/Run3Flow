import numpy as np
import ROOT
import pickle

import sys
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

data = {
        "Run2":ROOT.TFile("published/HEPData-ins1666817-v1-5TeVPbPb_flow.root","read"),
		"Run3":ROOT.TFile("data/fout_FT0C_v2_cos_v2.root","read"),
};

plotParams = {
        "Run2":{"color":"black","fmt":"s",'fillstyle':'none',"label":"5.02 TeV"},
		"Run3":{"color":"red","fmt":"o",'fillstyle':'none',"label":"5.36 TeV"}
        #"jyu_hydro":{"plotType":"theory","color":"orange","linestyle":"--","alpha":0.7,"label":"{T\\raisebox{-.5ex}{R}ENTo}+VISH(2+1)+UrQMD"},
};

cent =    [0,5,10,15,20,25,30,35,40,45,50,60,65,70,75,80];
icHEP =   [0,1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9, 9];

plot = JPyPlotRatio.JPyPlotRatio(panels=(3,3), #number of panels (each with ratio)
	xlabel="$p_\\mathrm{T}$", #create a 3x3 grid plot
	panelLabel={i: "{}--{} %".format(cent[i],cent[i+1]) for i in range(0,9)}, #label the panels v_n
	#panelScaling={1:2.0,2:3.0}, #add scaling to some of the panels (panel index : scale factor)
	systPatchWidth=0.08,
    ratioBounds ={0:(0.5,1.5),1:(0.1,2),2:(0.1,7)},
	ratioSystPlot=True,
	ylabel="$v_2$",disableRatio=[1,2]);

plots = {};
gr_vn = {};
for ic in range(0,9):
        #print("gr_v{}".format(i+2));
	if(ic<6):
		gr_vn[(0,ic)] = data["Run2"].Get("Table 2/Graph1D_y{}".format(icHEP[ic]))
		plots[0] = plot.Add(ic,gr_vn[(0,ic)],**plotParams["Run2"]);
	gr_vn[(1,ic)] = data["Run3"].Get("gPt_{}".format(ic));
	plots[1] = plot.Add(ic,gr_vn[(1,ic)],**plotParams["Run3"]);
	#plot.Ratio(plots["Run3"],plots["Run2"]);
	
	
	#plot.Ratio(hv2,hv2);

ax1 = plot.GetAxes(8); #get the last panel handle and add a long label manually
#ax1.text(0.05,0.32,"Pb-Pb $\\sqrt{s_\\mathrm{NN}}=5.02\\,\\mathrm{TeV}$\n$0.4<|\\eta|<0.8$, $0.2<p_\\mathrm{T}<5.0\\,\\mathrm{GeV/c}$",horizontalalignment="left",verticalalignment="center",transform=ax1.transAxes,size=9);

plot.Plot();
plot.Save("figs/vn_ptdiff_run3.pdf");
plot.Show();

