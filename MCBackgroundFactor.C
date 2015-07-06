
double	MCBackgroundFactor(const TFile* acquMcSignalFile, const TFile* acquMcBGFile, const TFile* out)
{
	TCanvas*	can = new TCanvas("canMCBackgroundFactor", "MCBackgroundFactor", 1500, 800);
    can->Divide(2, 1);
    
	TTree*	treeSignal	= (TTree*)acquMcSignalFile->Get("tagger");
	TTree*	treeBG		= (TTree*)acquMcBGFile->Get("tagger");
	
	cout << "calculating MCBackgroundFactor." << endl;
	cout << "acquMcSignalCount     acquMcBGCount     MCBackgroundFactor" << endl;
	
	can->cd(1);
	treeSignal->Draw("taggedChannel>>acquMcSignalCount");
	TH1D*	acquMcSignalCount = (TH1D*)gDirectory->Get("acquMcSignalCount");
	cout << acquMcSignalCount->GetBinContent(1) << "                  ";
	
	can->cd(2);
	treeBG->Draw("taggedChannel>>acquMcBGCount");
	TH1D*	acquMcBGCount = (TH1D*)gDirectory->Get("acquMcBGCount");
	cout << acquMcBGCount->GetBinContent(1) << "             ";
	
	out->cd();
	can->Write();
	
	double	res	= (3*acquMcSignalCount->GetBinContent(1))/(0.082*acquMcBGCount->GetBinContent(1));
	cout << res << endl;
	return res;
}

void	MCBackgroundFactor(const char* acquMcSignalFileName, const char* acquMcBGFileName)
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
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	MCBackgroundFactor(acquMcSignalFile, acquMcBGFile, out);
}


void	MCMCBackgroundFactor(const char* dir)
{
	std::strstream	acquMcSignalFileName;
	std::strstream	acquMcBGFileName;
	
	acquMcSignalFileName << dir << "/Acqu_g4_sim_etap.root";
	acquMcBGFileName << dir << "/Acqu_g4_sim_3pi0.root";
	
	MCMCBackgroundFactor(acquMcSignalFileName.str().c_str(), acquMcBGFileName.str().c_str());
}


void	MCMCBackgroundFactor()
{
	MCMCBackgroundFactor("goatTrees/Sim");
}
