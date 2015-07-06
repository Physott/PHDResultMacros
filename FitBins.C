#include "macros/NewBestFit.C"


double	TaggE[48] = {
1577.31,
1573.83,
1570.43,
1567.06,
1563.71,
1560.39,
1557.11,
1553.8,
1550.53,
1547.24,
1543.98,
1540.7,
1537.44,
1534.17,
1530.91,
1527.66,
1524.4,
1521.13,
1517.89,
1514.63,
1511.36,
1508.11,
1504.86,
1501.62,
1498.39,
1495.16,
1491.94,
1488.71,
1485.49,
1482.26,
1479.03,
1475.8,
1472.58,
1469.36,
1466.13,
1462.89,
1459.67,
1456.4,
1453.17,
1449.92,
1446.66,
1443.42,
1440.17,
1436.94,
1433.7,
1430.49,
1427.32,
1424.42};

void TaggEff(Double_t* array, Double_t* darray)
{
	array[0]	= 48.0/100.0;
	array[1]	= 47.0/100.0;
	array[2]	= 47.0/100.0;
	array[3]	= 49.0/100.0;
	array[4]	= 51.0/100.0;
	array[5]	= 53.0/100.0;
	array[6]	= 54.0/100.0;
	array[7]	= 55.0/100.0;
	array[8]	= 56.0/100.0;
	array[9]	= 56.0/100.0;
	array[10]	= 58.0/100.0;
	array[11]	= 57.0/100.0;
	array[12]	= 59.0/100.0;
	array[13]	= 59.0/100.0;
	array[14]	= 60.0/100.0;
	array[15]	= 58.0/100.0;
	array[16]	= 61.0/100.0;
	array[17]	= 60.0/100.0;
	array[18]	= 60.0/100.0;
	array[19]	= 61.0/100.0;
	array[20]	= 63.0/100.0;
	array[21]	= 62.0/100.0;
	array[22]	= 64.0/100.0;
	array[23]	= 64.0/100.0;
	array[24]	= 66.0/100.0;
	array[25]	= 62.0/100.0;
	array[26]	= 65.0/100.0;
	array[27]	= 62.0/100.0;
	array[28]	= 66.0/100.0;
	array[29]	= 63.0/100.0;
	array[30]	= 65.0/100.0;
	array[31]	= 63.0/100.0;
	array[32]	= 65.0/100.0;
	array[33]	= 65.0/100.0;
	array[34]	= 65.0/100.0;
	array[35]	= 63.0/100.0;
	array[36]	= 66.0/100.0;
	array[37]	= 62.0/100.0;
	array[38]	= 68.0/100.0;
	array[39]	= 65.0/100.0;
	array[40]	= 66.0/100.0;
	array[41]	= 65.0/100.0;
	array[42]	= 65.0/100.0;
	array[43]	= 64.0/100.0;
	array[44]	= 66.0/100.0;
	array[45]	= 62.0/100.0;
	array[46]	= 66.0/100.0;
	array[47]	= 0;
	for(int i=0; i<48; i++)
		array[i] += 0.02;
	for(int i=0; i<48; i++)
		darray[i]	= 2.0/100.0;
}

struct	Value
{
	Double_t	value;
	Double_t	error;
};
struct	FitValuesSignal
{
	Value	factor;
	Value	mean;
	Value	sigma;
	Double_t	count;
};
struct	FitValuesData
{
	FitValuesSignal	signal;
	FitValuesSignal	bg;
};
struct	FitValuesBG
{
	FitValuesSignal	gaus;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
	}pol2;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
	}pol3;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
		Value	a4;
	}pol4;
	struct
	{
		Value	a0;
		Value	a1;
		Value	a2;
		Value	a3;
		Value	a4;
		Value	a5;
	}pol5;
	Double_t	count;
};
struct	FitValues
{
	FitValuesSignal	signal;
	FitValuesBG		bg;
};



void	FitBinsSimGaus(const TH1D* data, int channel, int size, FitValues& fitValues, const bool signal)
{
	TF1*	fit = new TF1("fitfktSignal", "gaus(0)", 850, 1050);
	if(signal)
	{
		data->SetLineColor(kRed);
		fit->SetLineColor(kMagenta);
		fit->SetParameters(data->GetMaximum(), 957, 10);
		fit->SetParLimits(1, 950, 965);
		fit->SetParLimits(2, 5, 20);
	}
	else
	{
		data->SetLineColor(kGreen);
		fit->SetLineColor(kGreen);
		fit->SetParameters(data->GetMaximum(), 930, 50);
		fit->SetParLimits(1, 850, 1050);
		fit->SetParLimits(2, 25, 250);
	}
	fit->SetParLimits(0, 0, data->GetMaximum()*1.5);
	
	TString	title("E_{#gamma} (");
	if(channel==0)
	{
		title.Append(TString().Itoa(int(TaggE[size]-((TaggE[size]-TaggE[size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[0]+((TaggE[0]-TaggE[1])/2)), 10));
		title.Append(" MeV) ");
	}
	else if(channel+size>=47)
	{
		title.Append(TString().Itoa(int(TaggE[47]-((TaggE[46]-TaggE[47])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	else
	{
		title.Append(TString().Itoa(int(TaggE[channel+size]-((TaggE[channel+size]-TaggE[channel+size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	data->SetTitle(title.Data());
	data->GetXaxis()->SetTitle("IM 6#gamma [MeV]");	
	data->Fit(fit, "QR0");
	data->SetStats(0);		
    data->SetAxisRange(850, 1050);
	data->Draw("SAME");
	fit->Draw("SAME");
	
	if(signal)
	{
		fitValues.signal.factor.value	= fit->GetParameter(0);
		fitValues.signal.factor.error	= fit->GetParError(0);
		fitValues.signal.mean.value		= fit->GetParameter(1);
		fitValues.signal.mean.error		= fit->GetParError(1);
		fitValues.signal.sigma.value	= fit->GetParameter(2);
		fitValues.signal.sigma.error	= fit->GetParError(2);
		fitValues.signal.count			= data->Integral();
	}
	else
	{
		fitValues.bg.gaus.factor.value	= fit->GetParameter(0);
		fitValues.bg.gaus.factor.error	= fit->GetParError(0);
		fitValues.bg.gaus.mean.value	= fit->GetParameter(1);
		fitValues.bg.gaus.mean.error	= fit->GetParError(1);
		fitValues.bg.gaus.sigma.value	= fit->GetParameter(2);
		fitValues.bg.gaus.sigma.error	= fit->GetParError(2);
		fitValues.bg.count				= data->Integral();
	}
}
void	FitBinsSimPol(const TH1D* data, int channel, int size, FitValues& fitValues)
{
	TF1*	fit2 = new TF1("fitfktBG2", "pol2(0)", 850, 1050);
	TF1*	fit3 = new TF1("fitfktBG3", "pol3(0)", 850, 1050);
	TF1*	fit4 = new TF1("fitfktBG4", "pol4(0)", 850, 1050);
	TF1*	fit5 = new TF1("fitfktBG5", "pol5(0)", 850, 1050);
	fit2->SetLineColor(kGreen);
	fit3->SetLineColor(kOrange);
	fit4->SetLineColor(kYellow);
	fit5->SetLineColor(kCyan);
	fit2->SetParameters(0, 0, 0);
	fit3->SetParameters(0, 0, 0, 0);
	fit4->SetParameters(0, 0, 0, 0, 0);
	fit5->SetParameters(0, 0, 0, 0, 0, 0);
	
	data->Fit(fit2, "QR0");
	data->Fit(fit3, "QR0");
	data->Fit(fit4, "QR0");
	data->Fit(fit5, "QR0");
	fit2->Draw("SAME");
	fit3->Draw("SAME");
	fit4->Draw("SAME");
	fit5->Draw("SAME");
	
	fitValues.bg.pol2.a0.value	= fit2->GetParameter(0);
	fitValues.bg.pol2.a0.error	= fit2->GetParError(0);
	fitValues.bg.pol2.a1.value	= fit2->GetParameter(1);
	fitValues.bg.pol2.a1.error	= fit2->GetParError(1);
	fitValues.bg.pol2.a2.value	= fit2->GetParameter(2);
	fitValues.bg.pol2.a2.error	= fit2->GetParError(2);
	
	fitValues.bg.pol3.a0.value	= fit3->GetParameter(0);
	fitValues.bg.pol3.a0.error	= fit3->GetParError(0);
	fitValues.bg.pol3.a1.value	= fit3->GetParameter(1);
	fitValues.bg.pol3.a1.error	= fit3->GetParError(1);
	fitValues.bg.pol3.a2.value	= fit3->GetParameter(2);
	fitValues.bg.pol3.a2.error	= fit3->GetParError(2);
	fitValues.bg.pol3.a3.value	= fit3->GetParameter(3);
	fitValues.bg.pol3.a3.error	= fit3->GetParError(3);
	
	fitValues.bg.pol4.a0.value	= fit4->GetParameter(0);
	fitValues.bg.pol4.a0.error	= fit4->GetParError(0);
	fitValues.bg.pol4.a1.value	= fit4->GetParameter(1);
	fitValues.bg.pol4.a1.error	= fit4->GetParError(1);
	fitValues.bg.pol4.a2.value	= fit4->GetParameter(2);
	fitValues.bg.pol4.a2.error	= fit4->GetParError(2);
	fitValues.bg.pol4.a3.value	= fit4->GetParameter(3);
	fitValues.bg.pol4.a3.error	= fit4->GetParError(3);
	fitValues.bg.pol4.a4.value	= fit4->GetParameter(4);
	fitValues.bg.pol4.a4.error	= fit4->GetParError(4);
	
	fitValues.bg.pol5.a0.value	= fit5->GetParameter(0);
	fitValues.bg.pol5.a0.error	= fit5->GetParError(0);
	fitValues.bg.pol5.a1.value	= fit5->GetParameter(1);
	fitValues.bg.pol5.a1.error	= fit5->GetParError(1);
	fitValues.bg.pol5.a2.value	= fit5->GetParameter(2);
	fitValues.bg.pol5.a2.error	= fit5->GetParError(2);
	fitValues.bg.pol5.a3.value	= fit5->GetParameter(3);
	fitValues.bg.pol5.a3.error	= fit5->GetParError(3);
	fitValues.bg.pol5.a4.value	= fit5->GetParameter(4);
	fitValues.bg.pol5.a4.error	= fit5->GetParError(4);
	fitValues.bg.pol5.a5.value	= fit5->GetParameter(5);
	fitValues.bg.pol5.a5.error	= fit5->GetParError(5);
}
void	FitBinsSim(const TH1D* dataSignal, const TH1D* dataBG, int channel, int size, FitValues& fitValues)
{
	cout << "channel: " << channel << endl;
	FitBinsSimGaus(dataSignal, channel, size, fitValues, true);
	FitBinsSimGaus(dataBG, channel, size, fitValues, false);
	FitBinsSimPol(dataBG, channel, size, fitValues);
}

void	FitBinsSim(const TH2D* dataSignal, const TH2D* dataBG, const TFile* out, std::vector<int>& ch, std::vector<int>& size, const double signalScale, FitValues* fitValues)
{
	TCanvas*	can = new TCanvas("canFitBinsSim", "FitBinsSim", 1500, 800);
	std::cout << "can size: " << TMath::Ceil(TMath::Sqrt(ch.size())) << std::endl;
    can->Divide(TMath::Ceil(TMath::Sqrt(ch.size())), TMath::Ceil(TMath::Sqrt(ch.size())));
    
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		TH1D*	dataSignalProjX	= ((TH1D*)dataSignal->ProjectionX(TString("SigBin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1))->Clone();
		dataSignalProjX->Scale(signalScale);
		FitBinsSim( dataSignalProjX,
					(TH1D*)dataBG->ProjectionX(TString("BGBin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1),
					ch[i],
					size[i],
					fitValues[i]);
	}
	
	TCanvas*	canMain = new TCanvas("canFitBinsSimMain", "FitBinsSimMain", 1500, 800);
    canMain->Divide(3,3);
    
	canMain->cd(1);
	FitBinsSimGaus((TH1D*)dataSignal->ProjectionX(TString("SignalBinSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()], true);
	canMain->cd(2);
	FitBinsSimGaus((TH1D*)dataBG->ProjectionX(TString("BGBinGausSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()], false);
	canMain->cd(3);
	FitBinsSimGaus((TH1D*)dataBG->ProjectionX(TString("BGBinGausSimMain2").Data(), 1, 47), 0, 46, fitValues[ch.size()], false);
	FitBinsSimPol((TH1D*)dataBG->ProjectionX(TString("BGBinPolSimMain").Data(), 1, 47), 0, 46, fitValues[ch.size()]);
	
	double*	x = new double[ch.size()];
	double*	dx = new double[ch.size()];
	canMain->cd(4);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			x[i]	= 0;
			for(int k=0; k<size[i]; k++)
				x[i]	+=	TaggE[ch[i]+k];
			x[i]	/= size[i];
			dx[i]	= 0;
			
			y[i]	= fitValues[i].signal.count;
			dy[i]	= TMath::Sqrt(fitValues[i].signal.count);
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries Signal");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(5);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].bg.count;
			dy[i]	= TMath::Sqrt(fitValues[i].bg.count);
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries BG");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(6);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.count/fitValues[i].bg.count;
			dy[i]	= y[i]*TMath::Sqrt((1/fitValues[i].signal.count) + (1/fitValues[i].bg.count));
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("Entries relation (signal/bg)");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(7);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.factor.value;
			dy[i]	= fitValues[i].signal.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal factor");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(8);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.mean.value;
			dy[i]	= fitValues[i].signal.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal mean");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(9);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitValues[i].signal.sigma.value;
			dy[i]	= fitValues[i].signal.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("signal sigma");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	out->cd();
	canMain->Write();
	can->Write();
}


void	FitBins(const TH1D* data, int channel, int size, const FitValues& fitValues, FitValuesData& fitResult)
{
	TF1*	fit = new TF1("fitfkt", "gaus(0) + gaus(3)", 850, 1050);
	fit->SetParameters(data->GetMaximum(), 957, fitValues.signal.sigma.value, 0, fitValues.bg.gaus.mean.value, fitValues.bg.gaus.sigma.value);
	//fit->SetParameters(data->GetMaximum(), fitValues.signal.mean.value, fitValues.signal.sigma.value, 0, fitValues.bg.gaus.mean.value, fitValues.bg.gaus.sigma.value);
	fit->SetParLimits(0, 0, 1.5*data->GetMaximum());
	fit->SetParLimits(1, 0.99*fitValues.signal.mean.value, 1.01*fitValues.signal.mean.value);
	fit->SetParLimits(2, 0.6*fitValues.signal.sigma.value, 1.3*fitValues.signal.sigma.value);
	fit->SetParLimits(3, 0, data->GetMaximum());
	fit->SetParLimits(4, 0.9*fitValues.bg.gaus.mean.value, 1.1*fitValues.bg.gaus.mean.value);
	fit->SetParLimits(5, 0.9*fitValues.bg.gaus.sigma.value, 1.1*fitValues.bg.gaus.sigma.value);
	
	TString	title("E_{#gamma} (");
	if(channel==0)
	{
		title.Append(TString().Itoa(int(TaggE[size]-((TaggE[size]-TaggE[size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[0]+((TaggE[0]-TaggE[1])/2)), 10));
		title.Append(" MeV) ");
	}
	else if(channel+size>=47)
	{
		title.Append(TString().Itoa(int(TaggE[47]-((TaggE[46]-TaggE[47])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	else
	{
		title.Append(TString().Itoa(int(TaggE[channel+size]-((TaggE[channel+size]-TaggE[channel+size+1])/2)), 10));
		title.Append(" MeV - ");
		title.Append(TString().Itoa(int(TaggE[channel]+((TaggE[channel-1]-TaggE[channel])/2)), 10));
		title.Append(" MeV) ");
	}
	data->SetTitle(title.Data());
	data->GetXaxis()->SetTitle("IM #gamma#gamma [MeV]");	
	data->Fit(fit, "QR0");
	data->SetStats(0);		
    data->SetAxisRange(850, 1050);
	data->Draw("SAME");
	fit->Draw("SAME");
	
	TF1*	fitSignal = new TF1("fitfktSignal", "gaus(0)", 850, 1050);
	TF1*	fitBG = new TF1("fitfktBG", "gaus(0)", 850, 1050);
	fitSignal->SetLineColor(kRed);
	fitBG->SetLineColor(kGreen);
	fitSignal->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2));
	fitBG->SetParameters(fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5));
	fitSignal->Draw("SAME");
	fitBG->Draw("SAME");
	
	fitResult.signal.factor.value	= fit->GetParameter(0);
	cout << fitResult.signal.factor.value << endl;
	fitResult.signal.factor.error	= fit->GetParError(0);
	fitResult.signal.mean.value		= fit->GetParameter(1);
	fitResult.signal.mean.error		= fit->GetParError(1);
	fitResult.signal.sigma.value	= fit->GetParameter(2);
	fitResult.signal.sigma.error	= fit->GetParError(2);
	fitResult.bg.factor.value	= fit->GetParameter(3);
	fitResult.bg.factor.error	= fit->GetParError(3);
	fitResult.bg.mean.value		= fit->GetParameter(4);
	fitResult.bg.mean.error		= fit->GetParError(4);
	fitResult.bg.sigma.value	= fit->GetParameter(5);
	fitResult.bg.sigma.error	= fit->GetParError(5);
}
void	FitBins(const TH2D* data, const TFile* out, std::vector<int>& ch, std::vector<int>& size, const double signalScale, FitValues* fitValues, FitValuesData* fitResult)
{
	TCanvas*	can = new TCanvas("canFitBins", "FitBins", 1500, 800);
    can->Divide(TMath::Ceil(TMath::Sqrt(ch.size())), TMath::Ceil(TMath::Sqrt(ch.size())));
    
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		FitBins((TH1D*)data->ProjectionX(TString("Bin").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1),
				ch[i],
				size[i],
				fitValues[i],
				fitResult[i]);
	}
	
	
	TCanvas*	canMain = new TCanvas("canFitBinsMain", "FitBinsMain", 1500, 800);
    canMain->Divide(3,3);
    
	canMain->cd(1);
	FitBins((TH1D*)data->ProjectionX(TString("SignalBinMain").Data(), 1, 38), 0, 37, fitValues[ch.size()], fitResult[ch.size()]);
	
	double*	x = new double[ch.size()];
	double*	dx = new double[ch.size()];
	canMain->cd(4);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			x[i]	= 0;
			for(int k=0; k<size[i]; k++)
				x[i]	+=	TaggE[ch[i]+k];
			x[i]	/= size[i];
			dx[i]	= 0;
			
			y[i]	= fitResult[i].signal.factor.value;
			dy[i]	= fitResult[i].signal.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("factor Signal");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(5);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].signal.mean.value;
			dy[i]	= fitResult[i].signal.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("mean Signal");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(6);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].signal.sigma.value;
			dy[i]	= fitResult[i].signal.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("sigma Signal");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(7);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.factor.value;
			dy[i]	= fitResult[i].bg.factor.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("factor bg");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(8);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.mean.value;
			dy[i]	= fitResult[i].bg.mean.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("mean bg");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	canMain->cd(9);
	{
		double*	y = new double[ch.size()];
		double*	dy = new double[ch.size()];
		
		for(int i=0; i<ch.size(); i++)
		{
			y[i]	= fitResult[i].bg.sigma.value;
			dy[i]	= fitResult[i].bg.sigma.error;
		}
		
        TGraphErrors*	graph	= new TGraphErrors(ch.size(), x, y, dx, dy);
        graph->SetTitle("sigma bg");
        graph->GetXaxis()->SetTitle("photon energy E_{#gamma} [MeV]");
		graph->Draw();
	}
	
	out->cd();
	can->Write();
	canMain->Write();
}

void	ShowBoth(const TH2D* data, const TH2D* dataSignal, const TH2D* dataBG, const TFile* out, std::vector<int>& ch, std::vector<int>& size, const double signalScale, FitValuesData* fitResult)
{
	TCanvas*	can = new TCanvas("canFitBinsBoth", "FitBinsBoth", 1500, 800);
    can->Divide(TMath::Ceil(TMath::Sqrt(ch.size())), TMath::Ceil(TMath::Sqrt(ch.size())));
    
    FitValues*		fitVal	= new FitValues[ch.size()+1];
    FitValuesData*	fitRes	= new FitValuesData[ch.size()+1];
    
	for(int i=0; i<ch.size(); i++)
	{
		can->cd(i+1);
		TH1D*	dataSignalProjX	= ((TH1D*)dataSignal->ProjectionX(TString("SigBinBoth").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1))->Clone();
		dataSignalProjX->Scale(signalScale);
		FitBinsSimGaus(	dataSignalProjX,
						ch[i],
						size[i],
						fitVal[i],
						true);
		FitBinsSimGaus(	(TH1D*)dataBG->ProjectionX(TString("BGBinBoth").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1),
						ch[i],
						size[i],
						fitVal[i],
						false);
		TH1D*	hist	= (TH1D*)data->ProjectionX(TString("BinBoth").Append(TString().Itoa(i,10)).Data(), ch[i]+1, ch[i]+size[i]+1);
		cout << fitVal[i].signal.factor.value << endl;
		cout << fitResult[i].signal.factor.value << endl;
		cout << fitVal[i].signal.factor.value/fitResult[i].signal.factor.value << endl;
		hist->Scale(fitVal[i].signal.factor.value/fitResult[i].signal.factor.value);
		FitBins(hist,
				ch[i],
				size[i],
				fitVal[i],
				fitRes[i]);
	}
}

void ReconstructEff(const TFile* dataSignal, const TFile* acquSignalFile, const TFile* out, TH1D*& slice6, TH1D*& slice7, TH1D*& RecEffTot6, TH1D*& RecEffTot7)
{
	TTree*	tagger	= (TTree*)acquSignalFile->Get("tagger");
	TTree*	tracks	= (TTree*)acquSignalFile->Get("tracks");
	
    Int_t   nTagged;
	Int_t	taggedChannel[128];
	Int_t	nTracks;
	
	tagger->SetBranchAddress("nTagged", &nTagged);
	tagger->SetBranchAddress("taggedChannel", taggedChannel);
	tracks->SetBranchAddress("nTracks",&nTracks);
	
	TH1D*	count	= new TH1D("count", "count", 48, 0, 48);
	TH1D*	count6	= new TH1D("count6", "count6", 48, 0, 48);
	TH1D*	count7	= new TH1D("count7", "count7", 48, 0, 48);
	
	for(int i=0; i<tracks->GetEntriesFast(); i++)
	{
		tagger->GetEntry(i);
		tracks->GetEntry(i);
		
		for(int t=0; t<nTagged; t++)
		{
			count->Fill(taggedChannel[t]);
			if(nTracks==6)
				count6->Fill(taggedChannel[t]);
			else if(nTracks==7)
				count7->Fill(taggedChannel[t]);
		}
	}
	
	TCanvas*	can = new TCanvas("canReconstructionEfficiency", "ReconstructionEfficiency", 1500, 800);
    can->Divide(3, 4);
    
    can->cd(2);
    count6->Sumw2();
    count6->Scale(1);
    count6->SetTitle("all simulated #eta' events 6 Hits");
	count6->GetXaxis()->SetTitle("tagger channel");
	count6->SetStats(0);
    count6->Draw();
    can->cd(3);
    count7->Sumw2();
    count7->Scale(1);
    count7->SetTitle("all simulated #eta' events 7 Hits");
	count7->GetXaxis()->SetTitle("tagger channel");
	count7->SetStats(0);
    count7->Draw();
    
    can->cd(5);
    TH2D*	data	= (TH2D*)dataSignal->Get("WithoutProton/fit4/TaggerBinning/IM_Bins");
    slice6	= (TH1D*)((TH1D*)data->ProjectionY("SigRecEff6"))->Clone();
    slice6->SetTitle("accepted simulated #eta' events 6 Hits");
	slice6->GetXaxis()->SetTitle("tagger channel");
	slice6->SetStats(0);
	slice6->Draw();
	
    can->cd(6);
    data	= (TH2D*)dataSignal->Get("WithProton/fitProton6/TaggerBinning/IM_Bins");
    slice7	= (TH1D*)((TH1D*)data->ProjectionY("SigRecEff7"))->Clone();
    slice7->SetTitle("accepted simulated #eta' events 7 Hits");
	slice7->GetXaxis()->SetTitle("tagger channel");
	slice7->SetStats(0);
	slice7->Draw();
	
	can->cd(8);
    TH1D*	RecEff6	= slice6->Clone();
    RecEff6->Divide(count6);
    RecEff6->SetTitle("Rec Eff 6 Hits");
	RecEff6->GetXaxis()->SetTitle("tagger channel");
	RecEff6->SetStats(0);
    RecEff6->Draw();
    
    can->cd(9);
    TH1D*	RecEff7	= slice7->Clone();
    RecEff7->Divide(count7);
    RecEff7->SetTitle("Rec Eff 7 Hits");
	RecEff7->GetXaxis()->SetTitle("tagger channel");
	RecEff7->SetStats(0);
	RecEff7->Draw();
	
	can->cd(7);
    TTree*	tt	= (TTree*)acquSignalFile->Get("tracks");
    tt->Draw("nTracks");
	
    can->cd(10);
    count->Sumw2();
    count->Scale(1);
    count->SetTitle("all simulated #eta' events all Hits");
	count->GetXaxis()->SetTitle("tagger channel");
	count->SetStats(0);
    count->Draw();
    
    can->cd(11);
    RecEffTot6	= (TH1D*)slice6->Clone();
    RecEffTot6->Divide(count);
    RecEffTot6->SetTitle("Rec Eff total 6 Hits");
	RecEffTot6->GetXaxis()->SetTitle("tagger channel");
	RecEffTot6->SetStats(0);
    RecEffTot6->Draw();
    
    can->cd(12);
    RecEffTot7	= (TH1D*)slice7->Clone();
    RecEffTot7->Divide(count);
    RecEffTot7->SetTitle("Rec Eff total 7 Hits");
	RecEffTot7->GetXaxis()->SetTitle("tagger channel");
	RecEffTot7->SetStats(0);
	RecEffTot7->Draw();
	
	out->cd();
	can->Write();
}


void ScalerCorrection(const TFile* data, const TFile* out, TH1D*& scaler)
{
	TCanvas*	can = new TCanvas("canScaler", "Scaler", 1500, 800);
    can->Divide(1, 2);
    
	can->cd(1);
	TH1D*	scalerUC	= (TH1D*)data->Get("EPT_ScalerT");
	scalerUC->Draw();
	can->cd(2);
	scaler	= (TH1D*)data->Get("EPT_ScalerCorT");
	scaler->Draw();	
	
	out->cd();
	can->Write();
}

void	FitBins(const TFile* dataFile, const TFile* mcSinalFile, const TFile* mcBGFile, const TFile* acquSignalFile, const TFile* out, const Int_t addedChannels, const char* histName, const double signalScale)
{
	std::vector<int>	ch;
	std::vector<int>	size;
	while(addedChannels*ch.size()<47)
	{
		if(47-(addedChannels*ch.size())>=addedChannels)
			size.push_back(addedChannels);
		else
			size.push_back(47-(addedChannels*ch.size()));
		ch.push_back(addedChannels*ch.size());
		/*if(ch.back()==0)
		{
			std::cout << "taggedEnergyMin: " << TaggE[0]<< std::endl;
			std::cout << "taggedEnergyMin: " << TaggE[1]<< std::endl;
			std::cout << "taggedEnergyMin: " << (TaggE[1]-TaggE[0])/2 << std::endl;
			std::cout << "taggedEnergyMin: " << TaggE[0]-((TaggE[1]-TaggE[0])/2) << std::endl;
			taggedEnergyMin.puch_back(TaggE[0]-((TaggE[1]-TaggE[0])/2));
		}
		else
			taggedEnergyMin.puch_back(TaggE[ch.back()]-((TaggE[ch.back()]-TaggE[ch.back()-1])/2));
			
		if(47-(addedChannels*ch.size())>=addedChannels)
			taggedEnergyMax.puch_back(TaggE[ch.back()]+((TaggE[ch.back()+1]-TaggE[ch.back()])/2));
		else
			taggedEnergyMax.puch_back(TaggE[46]+((TaggE[46]-TaggE[45])/2));*/
			
		std::cout << "starting number: " << ch.back() << "     channel count: " << size.back() << std::endl;
	}
    
    FitValues*		fitValues	= new FitValues[ch.size()+1];
    FitValuesData*	fitResult	= new FitValuesData[ch.size()+1];
	TH2D*		dataSignal	= (TH2D*)mcSinalFile->Get(histName);
	TH2D*		dataBG		= (TH2D*)mcBGFile->Get(histName);
	//dataSignal->RebinX();
	dataBG->RebinX(10);
	dataBG->Scale(1.0/10.0);
    
    cout << endl << "FitBinsSim" << endl;
	FitBinsSim(dataSignal, dataBG, out, ch, size, signalScale, fitValues);
	
	TH2D*		data		= (TH2D*)dataFile->Get(histName);
	data->RebinX(5);
	data->Scale(1.0/5.0);
    cout << endl << "FitBins" << endl;
	FitBins(data, out, ch, size, signalScale, fitValues, fitResult);
	
	
    cout << endl << "ShowBoth" << endl;
	ShowBoth(data, dataSignal, dataBG, out, ch, size, signalScale, fitResult);
	
	
	TH1D*	simCount6;
	TH1D*	simCount7;
	TH1D*	recEff7;
	TH1D*	recEff6;
	TH1D*	recEff7;
	ReconstructEff(mcSinalFile, acquSignalFile, out, simCount6, simCount7, recEff6, recEff7);
	
	TH1D*	scaler;
	ScalerCorrection(dataFile, out, scaler);
	
	TCanvas*	can = new TCanvas("canResult7", "Result7", 1500, 800);
    can->Divide(2, 4);
	
	
	can->cd(2);
	TH1D*	count7	= new TH1D("resultCount7", "resultCount7", 48, 0, 48);
	count7->Sumw2();
	for(int i=0; i<ch.size(); i++)
	{
		count7->SetBinContent((addedChannels*i)+1, fitResult[i].signal.factor.value * fitResult[i].signal.sigma.value * TMath::Sqrt(TMath::Pi()));
		double	help1 = fitResult[i].signal.factor.error * fitResult[i].signal.sigma.value * TMath::Sqrt(TMath::Pi());
		double	help2 = fitResult[i].signal.factor.value * fitResult[i].signal.sigma.error * TMath::Sqrt(TMath::Pi());
		count7->SetBinError((addedChannels*i)+1, TMath::Sqrt((help1*help1)+(help2*help2)));
	}
	count7->SetLineColor(kBlack);
	count7->GetXaxis()->SetTitle("tagger channel");
	count7->SetTitle("count #eta'");
	count7->Draw();
	TH1D*	simCount7Rebin	= simCount7->Clone();
	//simCount7Rebin->Rebin(3, "width");
	//simCount7Rebin->GetXaxis()->Scale(1.0/3.0);
	simCount7Rebin->SetLineColor(kMagenta);
	simCount7Rebin->Scale(count7->GetMaximum()/simCount7Rebin->GetMaximum());
	simCount7Rebin->Draw("SAME");
	
	can->cd(3);
	TH1D*	recEff7rebin	= recEff7->Clone();
	recEff7rebin->Rebin(3, "width");
	recEff7rebin->SetTitle("reconstruction efficiency");
	recEff7rebin->GetXaxis()->SetTitle("tagger channel");
	recEff7rebin->Draw();
	
	can->cd(4);
	TH1D*	count7recEff	= count7->Clone();
	count7recEff->Divide(recEff7rebin);
	count7recEff->SetTitle("count #eta' reconstruct eff corrected");
	count7recEff->GetXaxis()->SetTitle("tagger channel");
	count7recEff->Draw();
	
	can->cd(5);
	double	tE[48];
	double	dtE[48];
	TaggEff(tE, dtE);
	TH1D*	hte	= new TH1D("TaggEFF", "TaggEFF", 48, 0, 48);
	for(int i=0; i<48; i++)
	{
		hte->SetBinContent(i+1, tE[i]);
		hte->SetBinError(i+1, dtE[i]);
	}
	hte->Rebin(3, "width");
	hte->SetTitle("tagging efficiency");
	hte->GetXaxis()->SetTitle("tagger channel");
	hte->GetXaxis()->SetRangeUser(0,47);
	hte->Draw();
	
	can->cd(6);
	TH1D*	count7recEfftaggEff	= count7recEff->Clone();
	count7recEfftaggEff->Divide(hte);
	count7recEfftaggEff->GetXaxis()->SetTitle("tagger channel");
	count7recEfftaggEff->SetTitle("count #eta' tagg eff corrected");
	count7recEfftaggEff->Draw();
	TH1D*	simCount7RebintaggEff	= simCount7Rebin->Clone();
	simCount7RebintaggEff->Scale(count7recEfftaggEff->GetMaximum()/simCount7RebintaggEff->GetMaximum());
	simCount7RebintaggEff->Draw("SAME");
	
	
	can->cd(7);
	scaler->SetTitle("scaler corrected");
	scaler->GetXaxis()->SetTitle("tagger channel");
	scaler->Draw();
	can->cd(8);
	TH1D*	result	= count7recEfftaggEff->Clone();
	result->Scale(1000000.0*4.2/0.08491);
	result->Scale(1.0/addedChannels);		//addedChannels
	result->Divide(scaler);
	result->SetTitle("cross section #eta'");
	result->GetXaxis()->SetTitle("tagger channel");
	result->GetYaxis()->SetTitle("cross-section #eta' [#mu b]");
	result->Draw();
	
	out->cd();
	can->Write();
	
	
	if(fitValues)		delete fitValues;
	if(fitResult)		delete fitResult;
}


void	FitBins(const char* dataFileName, const char* mcSignalFileName, const char* mcBGFileName, const char* acquSignalFileName, const Int_t addedChannels, const char* histName, const double	signalScale)
{
	TFile*	dataFile		= TFile::Open(dataFileName);
	if(!dataFile)
	{
		std::cout << "Can not open dataFile " << dataFileName << std::endl;
		return;
	}
	TFile*	mcSignalFile		= TFile::Open(mcSignalFileName);
	if(!mcSignalFile)
	{
		std::cout << "Can not open mcSignalFile " << mcSignalFileName << std::endl;
		return;
	}
	TFile*	mcBGFile			= TFile::Open(mcBGFileName);
	if(!mcBGFile)
	{
		std::cout << "Can not open mcBGFile " << mcBGFileName << std::endl;
		return;
	}
	TFile*	acquSignalFile			= TFile::Open(acquSignalFileName);
	if(!acquSignalFile)
	{
		std::cout << "Can not open scalerFile " << acquSignalFileName << std::endl;
		return;
	}
	TFile*	out				= TFile::Open("result.root", "RECREATE");
	if(!out)
	{
		std::cout << "Can not open output file result.root" << std::endl;
		return;
	}	
	
	FitBins(dataFile, mcSignalFile, mcBGFile, acquSignalFile, out, addedChannels, histName, signalScale);
}

void	FitBins(const char* dir, const Int_t addedChannels, const char* histName = "WithProton/fitProton6/TaggerBinning/IM_Bins")
{
	std::strstream	acquSignalFileName;
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
	double	signalScale	= BestFitSignalScale(mcSignalFile, mcBGFile);
	cout << "signalScale:   " << signalScale << endl;
	
	
	std::strstream	dataFileName;
	std::strstream	mcSignalFileName;
	std::strstream	mcBGFileName;
	std::strstream	scalerFileName;
	
	dataFileName << dir << "/Result_CB.root";
	mcSignalFileName << dir << "/Phys_g4_sim_etap_pi0pi0eta_00.root";
	mcBGFileName << dir << "/Phys_g4_sim_pi0pi0pi0_6g_00.root";
	
	cout << "FitBins:" << endl;
	FitBins(dataFileName.str().c_str(), mcSignalFileName.str().c_str(), mcBGFileName.str().c_str(), acquSignalFileName.str().c_str(), addedChannels, histName, signalScale);
}

void	FitBins()
{
	FitBins("goatTrees/Sim", 3, "WithProton/fitProton6/TaggerBinning/IM_Bins")
}
