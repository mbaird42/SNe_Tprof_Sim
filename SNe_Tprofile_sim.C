//==================================================//
//
// Supernova Time-Profile Simulator
// 
// code:     run_SNe_Tprofile_sim.C
// author:   Michael Baird
//           (m.baird@sussex.ac.uk, mbaird42@fnal.gov)
//
// Description: This macro runs the supernova time 
//   profile simulator.
//
// This code must be run in compiled mode!
// Example:
//          root -l run_SNe_Tprofile_sim.C+"(10.0,100,0.01,10,\"fake_SNe_input.root\",\"SNe_Tprofile_output.root\")"
//
//==================================================//

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom1.h"

#include <iostream>



// TODO:
//
// * put some double checks (for things like out of range values) into the findValueFrom2DHisto function
// * add some assert statements (about the size of the nu vectors etc...)
// * code in max number of tries for the findValue functions
// * remove the updateTP function and add "size of time window" to the trigger functions so that multiple time window sizes can be used in the same loop?



// ========== //
//
// ASSUMPTIONS:
//
// The following are a list of assumptions built into this code. The input histograms
// MUST adhere to these assumptions.
//
// 1. The SNe time-profile histogram assumes the number of expected events for a
//    distance of 10 kpc. This histo is renormalized by a factor of 1/r^2 where r is
//    the distance to the supernova (one of the macro inputs.)
//
// 2. The time profile histogram needs to have a linear binning for the X-axis.
//
// 3. The energy spectra (for both SNe and backgrounds) are interpreted as PDFs.
//    Therefore they must be normalized to have an area of 1.0.
//
// 4. The efficiency histograms MUST have the same X-axis binning as the energy
//    spectra, and they (obviously) must not exceed 1.0 for all bins.
//
// ========== //






// ================================================== //
//
// function to keep/reject an x-value from a 1D histo...
// (intended to be used for efficiencies)

bool keepValueFromEff(TH1F *h, double xValue, TRandom1 *rand) {

  // NOTE: It is assumed here that the efficiency is always a number between 0 and 1.0.

  bool keep = false;

  double xHeight = rand->Uniform(1.0);
  double content = h->GetBinContent(h->FindBin(xValue));
  if(xHeight <= content) {
    keep = true;
  }

  return keep;

}



// ================================================== //
//
// Function to pick out an x-value from a 1D histo.
// The 1D histo is assumed to be a normalized pdf.

double findValueFrom1DHisto(TH1F *h, TRandom1 *rand) {

  double value = 0.0;

  // extract info about the 1D histo:
  double xMin = h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst());
  double xMax = h->GetXaxis()->GetBinUpEdge (h->GetXaxis()->GetLast());
  double dX   = xMax-xMin;

  // find the max bin height:
  double maxHeightX = -1.0;
  for(int i = 1; i < h->GetNbinsX(); ++i) {
    double content = h->GetBinContent(i);
    if(content > maxHeightX) maxHeightX = content;
  }

  // pick a value:
  bool done = false;
  while(!done) {
    double xValue  = xMin + rand->Uniform(dX);
    double xHeight = rand->Uniform(maxHeightX);
    double content = h->GetBinContent(h->FindBin(xValue));
    if(xHeight <= content) {
      value = xValue;
      done  = true;
    }
  }

  return value;

}



// ================================================== //
//
// Function to pick out a y-value from a 2D histo.
// A vertical strip at X = xValue is taken to be a 1D pdf.

double findValueFrom2DHisto(TH2F *h, double xValue, TRandom1 *rand) {

  double value = 0.0;

  // extract info about the 2D histo:
  int    xBin = h->GetXaxis()->FindBin(xValue);
  double yMin = h->GetYaxis()->GetBinLowEdge(h->GetYaxis()->GetFirst());
  double yMax = h->GetYaxis()->GetBinUpEdge (h->GetYaxis()->GetLast());
  double dY   = yMax-yMin;

  // find the max bin height for the xValue column:
  double maxHeightY = -1.0;
  for(int j = 1; j < h->GetNbinsY(); ++j) {
    double content = h->GetBinContent(h->GetBin(xBin,j));
    if(content > maxHeightY) maxHeightY = content;
  }

  // pick a value:
  bool done = false;
  while(!done) {
    double yValue  = yMin + rand->Uniform(dY);
    double yHeight = rand->Uniform(maxHeightY);
    double content = h->GetBinContent(h->FindBin(xValue,yValue));
    if(yHeight <= content) {
      value = yValue;
      done  = true;
    }
  }

  return value;

}



// ================================================== //
//
// Function to clear out old entries from the TP lists.

void updateTPLists(std::vector<double> &times, std::vector<int> &nhits,
		   std::vector<int> &adc, std::vector<int> &truth,
		   double time, double winsize) {

  // NOTE: This function assumes that the elements in the list are
  //       in time order (which they should be unles the main code
  //       is drastically changed...)

  bool         useErase = false; // use the erase feature?
  unsigned int NErase   = 0;     // number of elements to erase
  for(unsigned int i = 0; i < times.size(); ++i) {
    if(times[i] < (time - winsize)) {
      useErase = true;
      NErase++;
    }
  }

  if(useErase) {
    times.erase(times.begin(),times.begin()+NErase);
    nhits.erase(nhits.begin(),nhits.begin()+NErase);
    adc  .erase(adc  .begin(),adc  .begin()+NErase);
    truth.erase(truth.begin(),truth.begin()+NErase);
  }

} 



// ================================================== //
//
// TRIGGER FUNCTION
//
// Simple trigger: just count # of TPs

bool TRIGGERcountTPs(std::vector<double> times, unsigned int Ntps) {

  bool trigger = false;

  if(times.size() >= Ntps) trigger = true;

  return trigger;

}



// ================================================== //
//
// TRIGGER FUNCTION
//
// Smart trigger: count # of TPs using NHit cuts

bool TRIGGERcountTPsWithNHitCut(std::vector<int> nhits, int loCut, int hiCut, unsigned int Ntps) {

  bool trigger = false;
  unsigned int count = 0;

  for(unsigned int i = 0; i < nhits.size(); i++) {
    if(nhits[i] >= loCut && nhits[i] <= hiCut) count++;
  }

  if(count >= Ntps) trigger = true;

  return trigger;

}






// ================================================== //
//
// main guts of the code...
//
// params:
//
//   * SNdist:      distance to supernova in kpc
//   * Nexpt:       number of experiments to generate
//   * WinSize:     size of sliding trigger window (in sec.)
//   * Ntps:        number of trigger primitives to trigger on
//   * inFileName:  string with input file name
//   * noiseScale:  by what factor to reduce/increase the noise

std::vector<double> generate_Nexpts(double SNdist, unsigned int Nexpt, double WinSize, int Ntps,
				    const char *inFileName, double noiseScale = 1.0) {

  // Print input info for the user:
  std::cout << "\n\n\n"
	    << "Running supernova time-profile simulator with input parameters:\n"
	    << "\tSNe distance:                      " << SNdist      << " kpc\n"
	    << "\tN experiments:                     " << Nexpt       << "\n"
	    << "\tsliding window size:               " << WinSize     << " sec\n"
	    << "\tnumber of TPs used for triggering: " << Ntps        << "\n"
	    << "\tinput file:                        " << inFileName  << "\n"
	    << "\tnoise scale factor:                " << noiseScale  << "\n"
	    << "\n\n";



  // open input file:
  TFile *inFile = new TFile(inFileName);



  // Get histograms from input file:
  TH1F *hTprofSN  = (TH1F*)inFile->FindObjectAny("hTprofSN");
  TH1F *hTprofBKG = (TH1F*)inFile->FindObjectAny("hTprofBKG");

  TH2F *hEspecSN  = (TH2F*)inFile->FindObjectAny("hEspecSN");
  TH1F *hEspecBKG = (TH1F*)inFile->FindObjectAny("hEspecBKG");

  TH2F *hNHitsSN  = (TH2F*)inFile->FindObjectAny("hNHitsSN");
  TH2F *hNHitsBKG = (TH2F*)inFile->FindObjectAny("hNHitsBKG");

  TH2F *hSumADCSN  = (TH2F*)inFile->FindObjectAny("hSumADCSN");
  TH2F *hSumADCBKG = (TH2F*)inFile->FindObjectAny("hSumADCBKG");

  TH1F *hEffSN  = (TH1F*)inFile->FindObjectAny("hEffSN");
  TH1F *hEffBKG = (TH1F*)inFile->FindObjectAny("hEffBKG");
  
  TH1F *hSplitTPprobSN  = (TH1F*)inFile->FindObjectAny("hSplitTPprobSN");
  TH1F *hSplitTPprobBKG = (TH1F*)inFile->FindObjectAny("hSplitTPprobBKG");
  
  
  
  // Extract info about the binning of the TprofSN histo. This info
  // is used to determine how to loop over time for each experiment.
  int    NtBins      = hTprofSN->GetNbinsX();
  double tBinLowEdge = hTprofSN->GetXaxis()->GetBinLowEdge(hTprofSN->GetXaxis()->GetFirst());
  double tBinUpEdge  = hTprofSN->GetXaxis()->GetBinUpEdge (hTprofSN->GetXaxis()->GetLast());

  // Renormalize the SN Tprof histo to the input distance:
  // (assumes the input histogram was scalled to 10.0 kpc!)
  hTprofSN->Scale(100.0/SNdist/SNdist);
  
  
  
  // Define a random number generator:
  TRandom1 *rand = new TRandom1();
  rand->SetSeed(0);



  // For handy printing of progress to the screen...
  double progress = 10.0;

  // Variables to keep track of true and false trigger rates.
  unsigned int NfalseTrigs   = 0;       // total number for all experiments
  unsigned int NtrueTrigs    = 0;       // total number for all experiments
  bool         falseTrig     = false;   // has a false trigger been issued?
  bool         trueTrig      = false;   // has a true trigger been issued?
  double       falseTrigTime = -10.0;   // time at which a false trigger has been issued

  // histos to be filled:



  //
  // The main loop:   Generate some fake SNe...
  //
  for(unsigned int n = 1; n <= Nexpt; ++n) {
    
    // Reset the trigger-tracker variables
    falseTrig     = false;
    trueTrig      = false;
    falseTrigTime = -10.0;

    // Print progress to the screen for impatient users...
    if((double)n/(double)Nexpt*100.0 >= progress) {
      std::cout << "\nlocal experiment set progress:\texperiment "
		<< n << " of " << Nexpt
		<< " (" << progress << "%)";
      progress += 10.0;
    }

    // Generate a random start time for the supernova between 1/4 and 1/2 way into the experiment.
    double SNstartTimeBin = (int)(NtBins/2 + rand->Uniform((int)(NtBins/2)));
    double SNstartTime    = (double)(SNstartTimeBin*(tBinUpEdge-tBinLowEdge)/NtBins);

    // background rate is assumed to be flat in time (thus the 1 bin histo)
    double NexpectedBKG = (hTprofBKG->GetBinContent(1))*noiseScale;



    // Define some vectors to keep running tallies of the TPs during the duration of the experiment
    std::vector<double> tpTimeList;
    std::vector<int>    tpNHitsList;
    std::vector<int>    tpSumADCList;
    std::vector<int>    tpTruthList;
    
    
    
    // Loop in time:  go in ticks of the hTprofSN binning for 2*(the size of that histo)
    for(int tBin = 1; tBin <= 2*NtBins; tBin++) {
      
      // Translate tBin (which is in "ticks") into a real time in seconds for ease of later calculations.
      double t = (double)(tBin*(tBinUpEdge-tBinLowEdge)/NtBins);

      // Pick a number of SNe nu events for this tick:
      unsigned int NevtsSN = 0;
      if(tBin >= SNstartTimeBin) {
	NevtsSN = rand->Poisson(hTprofSN->GetBinContent(1+tBin-SNstartTimeBin));
      }
      
      
      
      //
      // generate the info about the neutrinos (true E, NHits, Summed ADC)
      //
      
      // define vectors to hold the info about the individual neutrinos
      std::vector<double> nuE;
      std::vector<int>    nuNHits;
      std::vector<int>    nuSumADC;
      
      // loop over the number of SN nu events
      for(unsigned int i = 0; i < NevtsSN; ++i) {
	
	// pick the energy, number of hits, and the summed ADC for each true neutrino
	double energy = findValueFrom2DHisto(hEspecSN,t-SNstartTime,rand);
	int    nhits  = (int)findValueFrom2DHisto(hNHitsSN,energy,rand);
	int    sumadc = (int)findValueFrom2DHisto(hSumADCSN,energy,rand);
	
	bool keep = keepValueFromEff(hEffSN,energy,rand);

	// push the info into vectors if it passed the overall efficiency selector
	if(keep) {
	  nuE     .push_back(energy);
	  nuNHits .push_back(nhits);
	  nuSumADC.push_back(sumadc);
	}

      }
      
      
      
      // estimate pile up:
      // (loop over the list of true neutrinos, decide if some should be merged together, make new lists...)
      //
      // TODO:  This is not currently implemented!!!
      
      
      
      //
      // generate trigger primitives (TPs):
      //
      
      // define some vectors to store info about the TPs:
      std::vector<double> tpTime;   // the time of the TP
      std::vector<int>    tpNHits;  // the number of hits for this TP
      std::vector<int>    tpSumADC; // the summed ADC from all hits for this TP
      std::vector<int>    tpTruth;  // whether or not it was a true SN neutrino (use 1 for SN and 0 for background)
      
      // loop over the new list of true nu info, decide to make one or two TPs, make final vectors of TPs
      for(unsigned int i = 0; i < nuE.size(); ++i) {
	double prob = rand->Uniform(1.0);
	if(prob <= hSplitTPprobSN->GetBinContent(hSplitTPprobSN->FindBin(nuE[i]))) {
	  // split this nu into two TPs:
	  // split the NHits and ADC the same between the two TPs
	  double splitFrac = 0.1 + rand->Uniform(0.8); // don't let the split fraction for one be less than 10%
	  int nhit1 = (int)(nuNHits[i]*splitFrac);
	  int nhit2 = nuNHits[i] - nhit1;
	  int adc1  = (int)(nuSumADC[i]*splitFrac);
	  int adc2  = nuSumADC[i] - adc1;
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit1);
	  tpSumADC.push_back(adc1);
	  tpTruth .push_back(1);
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit2);
	  tpSumADC.push_back(adc2);
	  tpTruth .push_back(1);
	}
	else {
	  tpTime  .push_back(t);
	  tpNHits .push_back(nuNHits[i]);
	  tpSumADC.push_back(nuSumADC[i]);
	  tpTruth .push_back(1);
	}
      } // end loop over i
      
      
      
      
      
      // ========== //
      //
      // Repeat all of the above for the background events...
      //
      // ========== //
      
      // determine number of noise events in this time bin
      unsigned int NevtsBKG = 0;
      NevtsBKG = rand->Poisson(NexpectedBKG);
      
      
      //
      // generate the info about the background events (true E, NHits, Summed ADC)
      //
      
      // define vectors to hold the info about the individual neutrinos
      std::vector<double> bkgE;
      std::vector<int>    bkgNHits;
      std::vector<int>    bkgSumADC;
      
      // loop over the number of BKG events
      for(unsigned int i = 0; i < NevtsBKG; ++i) {
	
	// pick the energy, number of hits, and the summed ADC for each true neutrino
	double energy = findValueFrom1DHisto(hEspecBKG,rand);
	int   nhits  = (int)findValueFrom2DHisto(hNHitsBKG,energy,rand);
	int   sumadc = (int)findValueFrom2DHisto(hSumADCBKG,energy,rand);

	bool keep = keepValueFromEff(hEffBKG,energy,rand);
	
	// push the info into vectors if it passes the overall efficiency selector
	if(keep) {
	  bkgE     .push_back(energy);
	  bkgNHits .push_back(nhits);
	  bkgSumADC.push_back(sumadc);
	}
      }
      
      
      
      // estimate pile up:
      // (loop over the list of true neutrinos, decide if some should be merged together, make new lists...)
      //
      // TODO:  This is not currently implemented!!!
      
      
      
      //
      // generate trigger primitives (TPs):
      //
      
      // loop over the new list of true nu info, decide to make one or two TPs, make final vectors of TPs
      for(unsigned int i = 0; i < bkgE.size(); ++i) {
	double prob = rand->Uniform(1.0);
	if(prob <= hSplitTPprobBKG->GetBinContent(hSplitTPprobBKG->FindBin(bkgE[i]))) {
	  // split this nu into two TPs:
	  // split the NHits and ADC the same between the two TPs
	  double splitFrac = 0.1 + rand->Uniform(0.8); // don't let the split fraction for one be less than 10%
	  int nhit1 = (int)(bkgNHits[i]*splitFrac);
	  int nhit2 = bkgNHits[i] - nhit1;
	  int adc1  = (int)(bkgSumADC[i]*splitFrac);
	  int adc2  = bkgSumADC[i] - adc1;
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit1);
	  tpSumADC.push_back(adc1);
	  tpTruth .push_back(0);
	  
	  tpTime  .push_back(t);
	  tpNHits .push_back(nhit2);
	  tpSumADC.push_back(adc2);
	  tpTruth .push_back(0);
	}
	else {
	  tpTime  .push_back(t);
	  tpNHits .push_back(bkgNHits[i]);
	  tpSumADC.push_back(bkgSumADC[i]);
	  tpTruth .push_back(0);
	}
      } // end loop over i
      
      
      
      
      
      // Update the master list of TPs:
      //   * add new TPs to the list
      //   * wipe elements older than the look back time
      for(unsigned int i = 0; i < tpTime.size(); ++i) {
	tpTimeList  .push_back(tpTime[i]);
	tpNHitsList .push_back(tpNHits[i]);
	tpSumADCList.push_back(tpSumADC[i]);
	tpTruthList .push_back(tpTruth[i]);
      }
      updateTPLists(tpTimeList, tpNHitsList, tpSumADCList, tpTruthList, t, WinSize);
      
      
      
      //
      // Call the trigger function. For now, just call one trigger function at a time (comment in the one you want...)
      //
      bool trigger = false;

      // trigger 00:
      trigger = TRIGGERcountTPs(tpTimeList,Ntps);

      // trigger 01:
      // trigger = TRIGGERcountTPsWithNHitCut(tpNHitsList, 10, 50, Ntps);

      // If the trigger is after the true start of the SN && not more than 50 msec after the burst starts,
      // then consider the trigger successful.
      if(trigger) {
	    
	if((t-SNstartTime) >= 0.0 && (t-SNstartTime) <= 0.05) {
	  if(!trueTrig) {
	    /*
	    std::cout << "\n\n\nTRUE Trigger issued for experiment " << n << "\n"
		      << "at t = " << t << "\n"
		      << "SNe start time = " << SNstartTime
		      << "\n\n\n";
	    */
	    trueTrig = true;
	    NtrueTrigs++;
	  }
	} // end if(time is close to true burst time) block
	else {
	  // if a false trigger hasn't been issued in the last WinSize and we are in the first 1/4 of the experiment, then let a false trigger be issued
	  if(!falseTrig && t < 0.5*(tBinUpEdge-tBinLowEdge)) {
	    /*
	    std::cout << "\n\n\nFALSE Trigger issued for experiment " << n << "\n"
		      << "at t = " << t << "\n"
		      << "SNe start time = " << SNstartTime
		      << "\n\n\n";
	    */
	    falseTrig = true;
	    NfalseTrigs++;
	    falseTrigTime = t;
	  }
	} // end else block
	
      } // end if(trigger) block   

      // Reset the false trigger flag if one window size has elapsed since the last false trigger.
      if(t-falseTrigTime >= WinSize) {
	falseTrig = false;
      }
      
    } // end loop over time t
    
    // Fill some plots with the experimental results
    
    
  } // end loop over n (number of experiments)


  // Compute and return final numbers:
  double eff = (double)NtrueTrigs/(double)Nexpt;
  double ftr = (double)NfalseTrigs/((double)Nexpt*(0.5*(tBinUpEdge-tBinLowEdge)))*3.1557e7;

  std::cout << "\n\n"
	    << "\nEfficiency               = " << eff
	    << "\nnumber of false triggers = " << NfalseTrigs
	    << "\nfalse trigger rate       = " << ftr << " per year"
	    << "\n\n\n";

  std::vector<double> output;

  output.push_back(eff);
  output.push_back((double)NfalseTrigs);
  output.push_back(ftr);

  return output;

}



// ================================================== //
//
// wrapper to plot efficency vs. distance:
//

void eff_vs_dist(const char* inFileName, const char* outFileName,
		 unsigned int Nexpt, double WinSize, int Ntps, double noiseScale = 1.0) {

  // open the output file:
  TFile *outfile = TFile::Open(outFileName,"RECREATE");

  // book the efficiency histo:
  TH1F *hEff_vs_dist = new TH1F("hEff_vs_dist","SNe Trigger efficiency vs. SN distance;distance [kpc];Eff.",
				160,10.0,90.0);

  // loop over distances:
  for(int i = 1; i <= hEff_vs_dist->GetNbinsX(); ++i) {
    
    double SNdist = hEff_vs_dist->GetBinCenter(i);

    std::vector<double> output = generate_Nexpts(SNdist,Nexpt,WinSize,Ntps,inFileName,noiseScale);

    hEff_vs_dist->SetBinContent(i,output[0]);
    
    std::cout << "\n\n"
	      << "TOTAL progress:\t "
	      << 100.0*(double)i/((double)(hEff_vs_dist->GetNbinsX())) << "%"
	      << "\n\n";

  }
  
  // write output to disk:
  outfile->Write();

}



// ================================================== //
//
// wrapper to plot fake rate vs. noise rate
//

void fakeRate_vs_noiseRate(const char* inFileName, const char* outFileName,
			   double SNdist, unsigned int Nexpt, double WinSize, int Ntps) {

  // open the output file:
  TFile *outfile = TFile::Open(outFileName,"RECREATE");

  // book the fake rate histos:
  TH1F *hNFakeTrigs_vs_noiseRate = new TH1F("hNFakeTrigs_vs_noiseRate",
					    "Number of false triggers vs. noise/background scale;scale factor;number of triggers",
					    100,1.0,101.0);
  TH1F *hFakeRate_vs_noiseRate = new TH1F("hFakeRate_vs_noiseRate",
					  "Fake Trigger Rate vs. noise/background scale;scale factor;rate [triggers per year]",
					  100,1.0,101.0);

  // loop over distances:
  for(int i = 1; i <= hNFakeTrigs_vs_noiseRate->GetNbinsX(); ++i) {
    
    double noiseScale = hNFakeTrigs_vs_noiseRate->GetBinCenter(i);

    std::vector<double> output = generate_Nexpts(SNdist,Nexpt,WinSize,Ntps,inFileName,noiseScale);

    hNFakeTrigs_vs_noiseRate->SetBinContent(i,output[1]);
    hFakeRate_vs_noiseRate  ->SetBinContent(i,output[2]);

    std::cout << "\n\n"
	      << "TOTAL progress:\t "
	      << 100.0*(double)i/((double)(hNFakeTrigs_vs_noiseRate->GetNbinsX())) << "%"
	      << "\n\n";
    
  }
  
  // write output to disk:
  outfile->Write();

}



// ================================================== //
//
// help function:
//

void help() {

  std::cout << "\n\n"
	    << "===== HELP ====="
	    << "\n"
	    << "\nAvailable Functions:"
	    << "\n\n\n"

	    << "\n=========="
	    << "\nhelp()"
	    << "\n"
	    << "\n* Print this message."
	    << "\n* params: (none)"
	    << "\n\n\n"

	    << "\n=========="
	    << "\ngenerate_Nexpts(SNdist, Nexpt, WinSize, Ntps, \"inFileName\", noiseScale = 1.0)"
	    << "\n"
	    << "\n* Main loop to generate N experiments under specified conditions, returns a vector"
	    << "\n  of doubles (efficency, number of false triggers, fake trigger rate.) - Only"
	    << "\n  intended to be used by the wrapper functions listed below."
	    << "\n* params:"
	    << "\n   - SNdist     = distance to supernova in kpc"
	    << "\n   - Nexpt      = number of experiments to generate"
	    << "\n   - WinSize    = size of sliding trigger window in seconds"
	    << "\n   - Ntps       = number of trigger primitives to trigger on"
	    << "\n   - inFileName = file name with input histograms"
	    << "\n   - noiseScale = (optional, defaults to 1.0) adjust the assumed background rate"
	    << "\n\n\n"

	    << "\n=========="
	    << "\neff_vs_dist(\"inFileName\", \"outFileName\", Nexpt, WinSize, Ntps, noiseScale = 1.0)"
	    << "\n"
	    << "\n* Wrapper function for running generate_Nexpts with different distances. Generates"
	    << "\n  a histo of efficency vs. distance."
	    << "\n* params: (same as above plus...)"
	    << "\n   - outFileName = file name for writing output histogram(s)"
	    << "\n\n\n"

	    << "\n=========="
	    << "\nfakeRate_vs_noiseRate(\"inFileName\", \"outFileName\", SNdist, Nexpt, WinSize, Ntps)"
	    << "\n"
	    << "\n* Wrapper function for running generate_Nexpts with different noise rates. Generates"
	    << "\n  a histo of fake trigger rate and number of fake triggers vs. noise rate."
	    << "\n* params: (same as above)"
	    << "\n\n\n";

}

// ================================================== //
//
// primary function:
//
// Prints info on start up about what commands can be run.

void SNe_Tprofile_sim() {

  std::cout << "\n\n"
	    << "\n\t" << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
	    << "\n\t" << "*                                                                             *"
	    << "\n\t" << "*  Supernova Time Profile DAQ Simulator:                                      *"
	    << "\n\t" << "*  author:  Michael Baird (mbaird42@fnal.gov)                                 *"
	    << "\n\t" << "*                                                                             *"
	    << "\n\t" << "*  Welcome to the root command line interface for the supernova time profile  *"
	    << "\n\t" << "*  simulator! For a list of available commands, type 'help()'.                *"
	    << "\n\t" << "*                                                                             *"
	    << "\n\t" << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
	    << "\n\n\n\n";

}
