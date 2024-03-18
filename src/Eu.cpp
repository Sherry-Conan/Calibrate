#define Eu_cxx
#include "Eu.h"
#include <TStyle.h>
#include <TCanvas.h>
#include<iostream>
#include<sstream>
#include<algorithm>
#include <TSpectrum.h>
#include<vector>
#include<TGraph.h>
#include<TF1.h>
void Eu::GetBAndKmain(){
	Long64_t Nentry = fChain -> GetEntries(); 
	for(Int_t j = 0; j < 4; j++)	for(Int_t i = 0; i < 8; i++)	hGe[j][i] = new TH1I(Form("%02dh%02d", j, i), Form("%02dh%02d", j, i), 1500, 0, 6000); 
	for(Int_t i = 0; i < 4; i++)	hGe3[i] = new TH1I(Form("h%02d", i), Form("h%02d", i), 1500, 0, 6000); 
	for(Int_t i = 0; i < 10; i++)	hLaBr[i] = new TH1I(Form("hLaBr%02d", i), Form("hLaBr%02d", i), 1250, 0, 50000); 
	std::cout << Nentry << std::endl; 
	for(Long64_t Ientry = 0; Ientry < Nentry; Ientry++){
		if(Ientry % Int_t(1e7) == 0) std::cout<<Ientry<<std::endl; 
		fChain -> GetEntry(Ientry); 
		if(sid != 7 && cid == 0)
		{
			if(sid == 6) hGe3[ch - 4] -> Fill(evte); 
			else hGe[sid - 2][ch - 4] -> Fill(evte); 
		}
		else	if(cid == 0)	hLaBr[ch] -> Fill(evte); 
	}

	for(Int_t i = 0; i < 4; i++)
	{
		for(Int_t j = 0; j < 8; j++)
		{
			HuntPeak(hGe[i][j], HPGe, Short_t(sizeof(HPGe) / sizeof(HPGe[0])), 6000, 2, 0.05); 

			TF1* f = new TF1("f", "pol1", 0, 6000); 
			TGraph* g = new TGraph(Short_t(sizeof(HPGe) / sizeof(HPGe[0])), x2, HPGe); //数组必须为double

			g -> Fit(f, "SQ+", "NR", 0, 6000); 
			kbGe[1][i][j] = f -> GetParameter(0); 
			kbGe[0][i][j] = f -> GetParameter(1); 
			std::cout << "sid = " << i + 2 << "\t" << "ch = " << j + 4 << std::endl; 
		}
		HuntPeak(hGe3[i], HPGe, Short_t(sizeof(HPGe) / sizeof(HPGe[0])), 6000, 2, 0.05); 
		std::cout << "sid = " << 6 << "\t" << "ch = " << i + 4 << std::endl; 

		TF1* f = new TF1("f", "pol1", 0, 6000); 
		TGraph* g = new TGraph(Short_t(sizeof(HPGe) / sizeof(HPGe[0])), x2, HPGe); //数组必须为double

		g -> Fit(f, "SQ+", "NR", 0, 6000); 
		kbGe3[1][i] = f -> GetParameter(0); 
		kbGe3[0][i] = f -> GetParameter(1); 
	}
	for(Short_t i = 0; i < 10; i++)
	{
		HuntPeak(hLaBr[i], LaBr, Short_t(sizeof(LaBr) / sizeof(LaBr[0])), 200, 6, 0.01); 

		TF1* f = new TF1("f", "pol1", 0, 60000); 
		TGraph* g = new TGraph(sizeof(LaBr) / sizeof(LaBr[0]), x2, LaBr); //数组必须为double

		g -> Fit(f, "SQ+", "NR", 0, 60000); 
		kbLaBr[1][i] = f -> GetParameter(0); 
		kbLaBr[0][i] = f -> GetParameter(1); 
		std::cout << "sid = " << 7 << "\t" << "ch = " << i << std::endl; 
	}

	TFile* f = new TFile(Form("Calibrate%04d.root", 3092), "recreate"); 
	for(Int_t j = 0; j < 4; j++)
	{
		TDirectory* Folder1 = (TDirectory*) f -> mkdir(Form("sid%02d", j + 2)); 
		for(Int_t i = 0; i < 8; i++)	Folder1 -> WriteObject(hGe[j][i], Form("h%02d", i + 4)); 
	}
	TDirectory* Folder2 = (TDirectory*) f -> mkdir(Form("sid%02d", 6)); 
	for(Int_t i = 0; i < 4; i++)
	{
		Folder2 -> WriteObject(hGe3[i], Form("h%02d", i + 4)); 
	}
	TDirectory* Folder3 = (TDirectory*) f -> mkdir(Form("sid%02d", 7)); 
	for(Int_t i = 0; i < 10; i++)	Folder3 -> WriteObject(hLaBr[i], Form("hLaBr%02d", i)); 

	f -> Close(); 
}

void Eu::HuntPeak(TH1I* h, const Double_t* y, Short_t YLong, Short_t PeakHigh, Short_t Sigma = 2, Double_t pho = 0.05)
{
	std::cout << YLong << std::endl; 
	TSpectrum* PeakHunt = new TSpectrum(100);//寻到最多的峰，一般是从高到低。
	Int_t NumPeak = PeakHunt -> Search(h, Sigma, "", pho); //分别为寻峰的一维直方图，方差，可选参数不好用，阈值，越低容忍度越高，最高为1，最低为0。

	Double_t *XPeak, *YPeak; 
	XPeak = PeakHunt -> GetPositionX(); 
	YPeak = PeakHunt -> GetPositionY(); 

	UShort_t PeakWhere[100]; 
	Short_t NumP = 0; 

	for(Short_t i = 0; i < NumPeak; i++)
	{
		if(YPeak[i] > PeakHigh)
		{
			PeakWhere[NumP] = XPeak[i]; 
			std::cout << XPeak[i] << "\t"; 
			NumP++; 
		}
	}
	std::cout << std::endl; 
	
	std::sort(PeakWhere, PeakWhere + NumP); 

	ChiSquare = 1000000000; 
	Recursion(0, y, YLong, 1, PeakWhere, NumP); 

	for(Short_t k = 0; k < YLong; k++)
	{
		std::cout<<x2[k]<<"\t"; 
	}
	std::cout<<"ChiSquare = "<<ChiSquare<<std::endl; 
}
void Eu::Recursion(Short_t i, const Double_t* y, Short_t YLong, Short_t n, UShort_t* PeakWhere, Short_t NumP)
{
	if(n != YLong)	for(Short_t j = i; j < NumP - YLong + n; j++)
	{
		x1[n - 1] = PeakWhere[j]; 
		Recursion(j + 1, y, YLong, n + 1,  PeakWhere, NumP); 
	}
	else	for(Short_t j = i; j < NumP - YLong + n; j++)
	{
		x1[n - 1] = PeakWhere[j]; 
		
		TF1* f = new TF1("f", "pol1", 0, 60000); 
		TGraph* g = new TGraph(YLong, x1, y); //数组必须为double

		g -> Fit(f, "SQ+", "NR", 0, 60000); 

		if(ChiSquare > f -> GetChisquare())
		{
			ChiSquare = f -> GetChisquare(); 
			for(Short_t k = 0; k < YLong; k++)
			{
				x2[k] = x1[k]; 
			}
		}
	}
}

void Eu::WriteOut(std::ofstream &GeFile, std::ofstream &LaBrFile)
{
	for(Short_t i = 0; i < 4; i++)
	{
		for (Short_t j = 0; j < 8; j++)
		{
			GeFile << kbGe[0][i][j] << "\t" << kbGe[1][i][j] << std::endl; 
		}
	}
	for (Short_t i = 0; i < 4; i++)
	{
		GeFile << kbGe3[0][i] << "\t" << kbGe3[1][i] << std::endl; 
	}
	for (Short_t i = 0; i < 10; i++)
	{
		LaBrFile << kbLaBr[0][i] << "\t" << kbLaBr[1][i] << std::endl; 
	}
}