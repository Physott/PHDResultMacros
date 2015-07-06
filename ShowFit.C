#include "macros/MCBackgroundFactor.C"



void	ShowHistogram(const TH1D& data, const TH1D& mcSignal, const TH1D& mcBG)
{
	Double_t	max	= data.GetMaximum();
	
	TH1D*	both	= mcSignal.Clone();
	both->Add(&mcBG);
	Double_t	scale	= max/both->GetMaximum();
	both->Scale(scale);
	both->SetLineColor(kYellow);
	both->Draw();
	
	mcSignal.Scale(scale);
	mcSignal.SetLineColor(kMagenta);
	mcSignal.Draw("SAME");
	
	mcBG.Scale(scale);
	mcBG.SetLineColor(kGreen);
	mcBG.Draw("SAME");
	
	data.SetLineColor(kBlack);
	data.Draw("SAME");
}


void	OpenHistogram(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TString& str, const double signalScale)
{
	TH1D*	data		= (TH1D*)dataFile->Get(str);
	TH1D*	mcSignal	= (TH1D*)mcSignalFile->Get(str);
	mcSignal->Scale(signalScale);
	TH1D*	mcBG		= (TH1D*)mcBGFile->Get(str);
	ShowHistogram(*data, *mcSignal, *mcBG);
}

void	ShowFitCheckProton(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	TCanvas*		can	= new TCanvas("CanCheckProton", "CheckProton", 1500, 800);
	can->Divide(4,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_prAng", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/raw/_CheckProton_raw_copl", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_prAng", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutCoplanarity/_CheckProton_cutCopl_copl", signalScale);
	can->cd(3);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_prAng", signalScale);
	can->cd(7);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutProtonAngle/_CheckProtoncutPrAng_copl", signalScale);
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_prAng", signalScale);
	can->cd(8);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/CheckProton/cutBoth/_CheckProton_cutBoth_copl", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowFitRaw(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanRaw", "Raw", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_MM", signalScale);
	can->cd(3)->SetLogz();
	TH2D*	TOF	= (TH2D*)dataFile->Get("WithProton/Raw/_Raw_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/Raw/_Raw_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowFitSubIM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanSubIM", "SubIM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_MM", signalScale);
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/SubIM_Cut/_subIMCut_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/SubIM_Cut/_subIMCut_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}
void	ShowFitMM(const TFile* dataFile, const TFile* mcSignalFile, const TFile* mcBGFile, const TFile* out, const double signalScale)
{
	can	= new TCanvas("CanMM", "MM", 1500, 800);
	can->Divide(3,2);
	
	can->cd(1);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_IM", signalScale);
	can->cd(2);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_MM", signalScale);
	can->cd(3)->SetLogz();
	TOF	= (TH2D*)dataFile->Get("WithProton/MM_Cut/_TOF");
	TOF->Draw("colz");
	can->cd(4);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub0IM", signalScale);
	can->cd(5);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub1IM", signalScale);
	can->cd(6);
	OpenHistogram(dataFile, mcSignalFile, mcBGFile, "WithProton/MM_Cut/_MMCut_sub2IM", signalScale);
	
	out->cd();
	can->Write();
}

void	ShowFit(const TFile* fitDataFile, const TFile* fitMcSignalFile, const TFile* fitMcBGFile, const TFile* out, const double signalScale)
{	
	ShowFitCheckProton(fitDataFile, fitMcSignalFile, fitMcBGFile, out, signalScale);
	ShowFitRaw(fitDataFile, fitMcSignalFile, fitMcBGFile, out, signalScale);
	ShowFitSubIM(fitDataFile, fitMcSignalFile, fitMcBGFile, out, signalScale);
	ShowFitMM(fitDataFile, fitMcSignalFile, fitMcBGFile, out, signalScale);
}

void	ShowFit(const char* acquMcSignalFileName, const char* acquMcBGFileName, const char* fitDataFileName, const char* fitMcSignalFileName, const char* fitMcBGFileName)
{
	TFile*	acquMcSignalFile	= TFile::Open(acquMcSignalFileName);
	if(!acquMcSignalFile)
	{
		std::cout << "Can not open acquMcSignalFile " << acquMcSignalFileName << std::endl;
		return;
	}
	TFile*	acquMcBGFile		= TFile::Open(acquMcBGFileName);
	if(!acquMcBGFile)
	{
		std::cout << "Can not open acquMcBGFile " << acquMcBGFileName << std::endl;
		return;
	}
	TFile*	fitDataFile		= TFile::Open(fitDataFileName);
	if(!fitDataFile)
	{
		std::cout << "Can not open fitDataFile " << fitDataFileName << std::endl;
		return;
	}
	TFile*	fitMcSignalFile	= TFile::Open(fitMcSignalFileName);
	if(!fitMcSignalFile)
	{
		std::cout << "Can not open fitMcSignalFile " << fitMcSignalFileName << std::endl;
		return;
	}
	TFile*	fitMcBGFile		= TFile::Open(fitMcBGFileName);
	if(!fitMcBGFile)
	{
		std::cout << "Can not open fitMcBGFile " << fitMcBGFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	double signalScale	= 1/MCBackgroundFactor(acquMcSignalFile, acquMcBGFile, out);
	
	ShowFit(fitDataFile, fitMcSignalFile, fitMcBGFile, out, signalScale);
}

void	ShowFit(const char* dir)
{
	std::strstream	acquMcSignalFileName;
	std::strstream	acquMcBGFileName;
	std::strstream	fitDataFileName;
	std::strstream	fitMcSignalFileName;
	std::strstream	fitMcBGFileName;
	
	acquMcSignalFileName << dir << "/Acqu_g4_sim_etap.root";
	acquMcBGFileName << dir << "/Acqu_g4_sim_3pi0.root";
	fitDataFileName << dir << "/Fit_g4_sim_etap.root";
	fitMcSignalFileName << dir << "/Fit_g4_sim_etap.root";
	fitMcBGFileName << dir << "/Fit_g4_sim_3pi0.root";
	
	ShowFit(acquMcSignalFileName.str().c_str(), acquMcBGFileName.str().c_str(), fitDataFileName.str().c_str(), fitMcSignalFileName.str().c_str(), fitMcBGFileName.str().c_str());
	
	/*std::strstream	acquSignalFileName;
	std::strstream	acquBGFileName;
	
	acquSignalFileName << dir << "/Acqu_g4_sim_etap_pi0pi0eta_00.root";
	acquBGFileName << dir << "/Acqu_g4_sim_pi0pi0pi0_6g_00.root";
	TFile*	mcSignalFile	= TFile::Open(acquSignalFileName.str().c_str());
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << acquSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile		= TFile::Open(acquBGFileName.str().c_str());
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << acquBGFileName << std::endl;
		return;
	}
	double	signalScale	= SignalScale(mcSignalFile, mcBGFile);
	cout << "signalScale:   " << signalScale << endl;
	
	
	std::strstream	dataFileName;
	std::strstream	mcSignalFileName;
	std::strstream	mcBGFileName;
	
	dataFileName << dir << "/Result_CB.root";
	mcSignalFileName << dir << "/Phys_g4_sim_etap_pi0pi0eta_00.root";
	mcBGFileName << dir << "/Phys_g4_sim_pi0pi0pi0_6g_00.root";
	
	cout << "FitBins:" << endl;
	ShowFit(dataFileName.str().c_str(), mcSignalFileName.str().c_str(), mcBGFileName.str().c_str(), signalScale);*/
}

void	ShowFit()
{
	ShowFit("goatTrees/Sim");
}
