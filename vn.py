import numpy as np
import ROOT
import pickle

import sys
sys.path.append("JPyPlotRatio");
import JPyPlotRatio

data = {
        "Run2":ROOT.TFile("published/lhc15o.root","read"),
		"Run3":ROOT.TFile("data/vn_run3.root","read"),
		"Run3_v2":ROOT.TFile("data/vn_run3_v2.root","read"),
		"Run3_v3":ROOT.TFile("data/vn_run3_pass4.root","read"),
		"Run3_v4":ROOT.TFile("data/vn_run3_pass4-JP.root","read"),
};

plotParams = {
        "Run2":{"color":"black","fmt":"s",'fillstyle':'none',"label":"5.02 TeV"},
		"Run3":{"color":"red","fmt":"o",'fillstyle':'none',"label":"5.36 TeV"},
        "Run3_v2":{"color":"red","fmt":"D",'fillstyle':'none',"label":"5.36 v2 TeV"},
        "Run3_v3":{"color":"blue","fmt":"*",'fillstyle':'none',"label":"5.36 pass4 TeV"},
        "Run3_v4":{"color":"cyan","fmt":"+",'fillstyle':'none',"label":"5.36 pass4 JP TeV"},
        #"jyu_hydro":{"plotType":"theory","color":"orange","linestyle":"--","alpha":0.7,"label":"{T\\raisebox{-.5ex}{R}ENTo}+VISH(2+1)+UrQMD"},
};


plot = JPyPlotRatio.JPyPlotRatio(panels=(3,3), #number of panels (each with ratio)
	xlabel="Centrality (%)", #create a 3x3 grid plot
	panelLabel={i: "$v_{{{}}}$".format(i+2) for i in range(0,8)}, #label the panels v_n
	panelScaling={1:2.0,2:3.0}, #add scaling to some of the panels (panel index : scale factor)
	systPatchWidth=0.08,
    ratioBounds ={0:(0.5,1.5),1:(0.1,2),2:(0.1,7)},
	ratioSystPlot=True,
	ylabel="$v_n$",disableRatio=[1,2]);

plots = {};
gr_vn = {};
for i in range(0,8):
	for s in list(data)[::-1]:
                #print("gr_v{}".format(i+2));
                gr_vn[(s,i)] = data[s].Get(("gr_v{}").format(i+2));
                plots[s] = plot.Add(i,gr_vn[(s,i)],**plotParams[s]);
	plot.Ratio(plots["Run3"],plots["Run2"]);
	plot.Ratio(plots["Run3_v2"],plots["Run2"]);
	plot.Ratio(plots["Run3_v3"],plots["Run2"]);
	
	#plot.Ratio(hv2,hv2);

ax1 = plot.GetAxes(8); #get the last panel handle and add a long label manually
ax1.text(0.05,0.32,"Pb-Pb $\\sqrt{s_\\mathrm{NN}}=5.02\\,\\mathrm{TeV}$\n$0.4<|\\eta|<0.8$, $0.2<p_\\mathrm{T}<5.0\\,\\mathrm{GeV/c}$",horizontalalignment="left",verticalalignment="center",transform=ax1.transAxes,size=9);

plot.Plot();
plot.Save("figs/vn.pdf");
plot.Show();

