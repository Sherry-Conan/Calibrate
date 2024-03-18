#ifndef Eu_h
#define Eu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <TH1I.h>
// Header file for the classes stored in the TTree if any.

const Double_t	HPGe[] = {244.5, 344.5, 356.5, 778.5}; //确定峰位，本实验主要集中在1000keV以下，故不对更高能量进行刻度，建议使用五个或以上，但不要超过10个
const Double_t	LaBr[] = {80.5, 121.5, 244.5, 352.5, 511, 778.5}; 
const Short_t yLong = std::max(sizeof(HPGe) / sizeof(HPGe[0]), sizeof(LaBr) / sizeof(LaBr[0])); 

class Eu {
public :
	TTree		*fChain = NULL;

	UShort_t	evte; 
	Short_t		ch; 
	Short_t		sid; 
	Short_t		cid; 

	Double_t	kbGe[2][4][8]; 
	TH1I*		hGe[4][8]; 
	TH1I*		hLaBr[10]; 
	Double_t	kbGe3[2][4]; 
	TH1I*		hGe3[4]; 
	Double_t	kbLaBr[2][10]; 

	Double_t		x1[yLong] = {0}; 
	Double_t		x2[yLong] = {0}; 
	Double_t		ChiSquare = 1000000000; 

   // List of branches
	TBranch		*b_evte;   //!
	TBranch		*b_ch;   //!
	TBranch		*b_sid;   //!
	TBranch		*b_cid;   //!

	Eu(Int_t FileNum); 

	void	Init(); 
	void	GetBAndKmain(); 
	void	WriteOut(std::ofstream &kbFile, std::ofstream &LaBrFile); 
	void	HuntPeak(TH1I* h, const Double_t* y, Short_t YLong, Short_t PeakHigh, Short_t Sigma, Double_t pho); 
	void	Recursion(Short_t i, const Double_t* y, Short_t YLong, Short_t n, UShort_t* PeakWhere, Short_t NumP); 
};

#endif

#ifdef Eu_cxx
Eu::Eu(Int_t FileNum){
	TFile *f = new TFile(Form("../../Data/rootfile/data_C1_3092_wave.root")); 
	fChain = (TTree*)f -> Get("tree"); 
	Init();
}

void Eu::Init(){
	fChain->SetBranchAddress("evte", &evte, &b_evte); 
	fChain->SetBranchAddress("ch", &ch, &b_ch); 
	fChain->SetBranchAddress("sid", &sid, &b_sid); 
	fChain->SetBranchAddress("cid", &cid, &b_cid); 
}
#endif // #ifdef Eu_cxx