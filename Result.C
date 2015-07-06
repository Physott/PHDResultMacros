#include "macros/ShowFit.C"



void	Result(const TFile* acquMcSignalFile, const TFile* acquMcBGFile, const TFile* fitDataFile, const TFile* fitMcSignalFile, const TFile* fitMcBGFile, const TFile* out)
{
	double backgroundFactor	= MCBackgroundFactor(acquMcSignalFile, acquMcBGFile, out);
	ShowFit(fitDataFile, fitMcSignalFile, fitMcBGFile, out, 1/backgroundFactor);
}

void	Result(const char* acquMcSignalFileName, const char* acquMcBGFileName, const char* fitDataFileName, const char* fitMcSignalFileName, const char* fitMcBGFileName)
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
	
	Result(acquMcSignalFile, acquMcBGFile, fitDataFile, fitMcSignalFile, fitMcBGFile, out);
}


void	Result(const char* dir)
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
	
	Result(acquMcSignalFileName.str().c_str(), acquMcBGFileName.str().c_str(), fitDataFileName.str().c_str(), fitMcSignalFileName.str().c_str(), fitMcBGFileName.str().c_str());
}


void	Result()
{
	Result("goatTrees/Sim");
}
