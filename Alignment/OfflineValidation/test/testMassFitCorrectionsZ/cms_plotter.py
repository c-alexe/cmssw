import ROOT as r
import cmsstyle as CMS
from ROOT import gStyle, TCanvas, TH1D, TLegend, kRed, kBlue

gStyle.SetOptStat(0)

'''
##############
# 2D HISTOGRAM
##############

# File reading 
f = r.TFile.Open('/eos/user/c/calexe/massfit_PostVFP_Iter0.root')
hist2d = f.Get("hcorin_0")

# Plotting
CMS.SetExtraText("")
CMS.SetCmsText("Private work (CMS data)")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75*0.76)
CMS.SetLumi("")

canv_2D = CMS.cmsCanvas('hist2d_root',0,72,0,72,'Internal A, #varepsilon, M (#eta)','Internal A, #varepsilon, M (#eta)',square=True,extraSpace=0.01,iPos=0,with_z_axis=True)

hist2d.Draw("same colz")

hdf = CMS.GetcmsCanvasHist(canv_2D)
hdf.GetXaxis().SetMaxDigits(2)

hdf.GetXaxis().SetTitleOffset(1.35)
hdf.GetXaxis().SetTitleSize(0.035)
hdf.GetXaxis().SetLabelSize(0.035)
hdf.GetXaxis().SetLabelOffset(0.005)
hdf.GetXaxis().SetMaxDigits(0)

hdf.GetYaxis().SetTitleOffset(1.35)
hdf.GetYaxis().SetTitleSize(0.035)
hdf.GetYaxis().SetLabelSize(0.035)
hdf.GetYaxis().SetLabelOffset(0.005)

#hist2d.GetZaxis().SetTitle("Events")
hist2d.GetZaxis().SetTitleOffset(1.35)
hist2d.GetZaxis().SetTitleSize(0.035)
hist2d.GetZaxis().SetLabelSize(0.035)
hist2d.GetZaxis().SetLabelOffset(0.005)

# Set the CMS official palette
#CMS.SetCMSPalette("")

# Allow to adjust palette position
#CMS.UpdatePalettePosition(hist2d, canv)

canv_2D.Print("corin_2016.pdf]")
'''

##############
# 1D HISTOGRAM
##############

# File reading 
f = r.TFile.Open('/eos/user/c/calexe/massscales_PostVFP_Iter0.root')
hist1d = f.Get("h_scales")

# Plotting
CMS.SetExtraText("")
CMS.SetCmsText("Private work (CMS data)")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75*0.76)
CMS.SetLumi("")

canv_1D = CMS.cmsCanvas('hist1d_root',0,21000,0,1.5,'bin no.','scale',square=False,extraSpace=0.01,iPos=0)

#canv_1D.SetEnergy(13)
#canv_1D.SetFrameBorderSize(1)

hist1d.Draw("same")

hdf1 = CMS.GetcmsCanvasHist(canv_1D)

hdf1.GetXaxis().SetTitleOffset(1.35)
hdf1.GetXaxis().SetTitleSize(0.035)
hdf1.GetXaxis().SetLabelSize(0.035)
hdf1.GetXaxis().SetLabelOffset(0.005)

hdf1.GetYaxis().SetTitleOffset(1.35)
hdf1.GetYaxis().SetTitleSize(0.035)
hdf1.GetYaxis().SetLabelSize(0.035)
hdf1.GetYaxis().SetLabelOffset(0.005)

canv_1D.Print("test1dnew.pdf]")

'''
####################
# 1D STACK HISTOGRAM
####################

# File reading
f = r.TFile.Open('/eos/user/c/calexe/massscales_SmearRealistic_toy0_Iter0_plots.root')

hist_d = f.Get("/postfit/h_data_10327")
hist_mc = f.Get("/postfit/h_prefit_10327")

# Plotting
CMS.SetExtraText("")
CMS.SetCmsText("Private work (CMS data)")
CMS.SetCmsTextFont(52)
CMS.SetCmsTextSize(0.75*0.76)
CMS.SetLumi("")

canv_1D = CMS.cmsCanvas('hist1d_root',70,110,0,700,'mll [GeV]','Events',square=True,extraSpace=0.01)

#canv_1D.SetEnergy(13)
#canv_1D.SetFrameBorderSize(1)

hist_d.SetLineColor(kBlue)
hist_mc.SetLineColor(kRed)
hist_d.Draw("hist same")
hist_mc.Draw("hist same")

hdf1 = CMS.GetcmsCanvasHist(canv_1D)

hdf1.GetXaxis().SetTitleOffset(1.35)
hdf1.GetXaxis().SetTitleSize(0.035)
hdf1.GetXaxis().SetLabelSize(0.035)
hdf1.GetXaxis().SetLabelOffset(0.005)

#hdf1.GetYaxis().SetMaximum(1000)
hdf1.GetYaxis().SetTitleOffset(1.35)
hdf1.GetYaxis().SetTitleSize(0.035)
hdf1.GetYaxis().SetLabelSize(0.035)
hdf1.GetYaxis().SetLabelOffset(0.005)

canv_1D.Print("mfit.pdf]")
'''