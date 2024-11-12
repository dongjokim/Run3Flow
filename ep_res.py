import numpy as np
import ROOT
import pickle

import sys
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

data = {
        "R2":ROOT.TFile("EPdata/EventPlaneResolution_psi2.root","read"),
		"R3":ROOT.TFile("EPdata/EventPlaneResolution_psi3.root","read"),
		"R4":ROOT.TFile("EPdata/EventPlaneResolution_psi4.root","read")
};

plotParams = {
        "R2":{"color":"black","fmt":"s",'fillstyle':'none',"label":"n = 2"},
		"R3":{"color":"red","fmt":"o",'fillstyle':'none',"label":"n = 3"},
        "R4":{"color":"blue","fmt":"D",'fillstyle':'none',"label":"n = 4"},
        #"jyu_hydro":{"plotType":"theory","color":"orange","linestyle":"--","alpha":0.7,"label":"{T\\raisebox{-.5ex}{R}ENTo}+VISH(2+1)+UrQMD"},
};
ylimits = [(-0.05,0.99)];

plot = JPyPlotRatio.JPyPlotRatio(panels=(1,1), #number of panels (each with ratio)
	xlabel="Centrality (%)", #create a 3x3 grid plot
	rowBounds=ylimits,  # for nrow
	disableRatio=[0],
	#panelLabel={i: "$v_{{{}}}$".format(i+2) for i in range(0,8)}, #label the panels v_n
	panelScaling={1:2.0,2:3.0}, #add scaling to some of the panels (panel index : scale factor)
	systPatchWidth=0.08,
    ratioBounds ={0:(0.5,1.5),1:(0.1,2),2:(0.1,7)},
	ratioSystPlot=True,
	ylabel="$R_n$",
	legendLoc=(0.81,0.83),legendSize=10,
	#disableRatio=[1,2]
	);

plots = {};
gr_vn = {};
for i in range(0,1):
	for s in list(data)[::-1]:
		#print("gr_v{}".format(i+2));
		gr_vn[(s,i)] = data[s].Get(("gRes"))
		
		plots[s] = plot.Add(i,gr_vn[(s,i)],**plotParams[s])
		
	
#plot.Ratio(plots["Run3"],plots["Run2"]);

ax1 = plot.GetAxes(0); #get the last panel handle and add a long label manually
ax1.text(0.05,0.92,"Pb-Pb $\\sqrt{s_\\mathrm{NN}}=5.36TeV$",horizontalalignment="left",verticalalignment="center",transform=ax1.transAxes,size=9);

plot.Plot();
plot.Save("figs/EP_res_Run3.pdf");
plot.Show();

