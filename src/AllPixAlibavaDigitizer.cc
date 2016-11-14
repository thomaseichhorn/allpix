/**
 *  Author:
 *    Thomas Eichhorn <thomas.eichhorn@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixAlibavaDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandGauss.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"

AllPixAlibavaDigitizer::AllPixAlibavaDigitizer(G4String modName, G4String hitsColName, G4String digitColName) : AllPixDigitizerInterface (modName)
{

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// threshold
	m_digitIn.thl = 0.;

	InitVariables();
}

void AllPixAlibavaDigitizer::InitVariables()
{
	// %/100 mean of charge to be shared
	sharemean = 0.225;

	// %/100 sigma of charge to be shared
	sharesig = 0.05;

	// 1- x %/100 by how much signals are scaled down if they are not in the centre of the pixel
	// 0.14132/44.479 1/um from fit in run000022 epi 100p unirr data
	// 0.02712/28.0768 1/um from fit in run000225 epi 100p 1.3e16 data
	distancereduce = 0.14132/44.479;

	// gain factor for signals going from electrons to adcs
	gain = 1/183.0;

	// 1/183 from electrons to adcs
	// was ca. 400 * [signal in MeV]

	// threshold for digitizer output
	threshold = 0.0;

	//
	// vars for output plots, not relevant for the simulation
	//
	
	// seedcut for plotting
	seedcut = 5.0;

	// clustercut for plotting
	clustercut = 2.5;

	// expected noise mean
	noisemean = 0.0;

	// expected noise sigma
	noisesigma = 5.0;

	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	NPixelX = gD->GetNPixelsX();
	NPixelY = gD->GetNPixelsY();
	pitchX=gD->GetPixelX()/um;
	pitchY=gD->GetPixelY()/um;

	rootoutputfile = new TFile("someoutputfile.root","RECREATE");

	histochargeleft = new TH1D("Output Charge Left Channel","Output Charge Left Channel",211,-10,200);
	histochargeleft->SetTitle("Output Charge Left Channel;ADCs;Entries");
	
	histochargemiddle = new TH1D("Output Charge Middle Channel","Output Charge Middle Channel",211,-10,200);
	histochargemiddle->SetTitle("Output Charge Middle Channel;ADCs;Entries");

	histochargeright = new TH1D("Output Charge Right Channel","Output Charge Right Channel",211,-10,200);
	histochargeright->SetTitle("Output Charge Right Channel;ADCs;Entries");
	
	histohitcount = new TH1I("Hit Count per Track","Hit Count per Track",11,-0.5,10.5);
	histohitcount->SetTitle("Hit Count per Track;Hits per Track;Entries");

	histohitcountoverthresh = new TH1I("Hit Count per Track over Digitizer Threshold","Hit Count per Track over Digitizer Threshold",11,-0.5,10.5);
	histohitcountoverthresh->SetTitle("Hit Count per Track over Digitizer Threshold;Hits per Track;Entries");

	histoclustersize = new TH1I("CShareclustersize","CShareclustersize",5,-0.5,4.5);
	histoclustersize->SetTitle("Expected Cluster Size from Charge Sharing;Cluster Size;Entries");

	histoloss_position = new TH1D("Charge Lost due to Track Position","Charge Lost due to Track Position",250,0,250);
	histoloss_position->SetTitle("Charge Lost due to Track Position;ADCs;Entries");

	histoloss_share = new TH1D("Charge Lost due to Charge Sharing","Charge Lost due to Charge Sharing",250,0,250);
	histoloss_share->SetTitle("Charge Lost due to Charge Sharing;ADCs;Entries");

	histoadcout = new TH1D("Track Induced Signal in ADCs","Track Induced Signal in ADCs",250,0,250);
	histoadcout->SetTitle("Track Induced Signal in ADCs;ADCs;Entries");

	histoeleout = new TH1D("Track Induced Signal in Electrons","Track Induced Signal in Electrons",250,0,40000);
	histoeleout->SetTitle("Track Induced Signal in Electrons;Electrons;Entries");

}

AllPixAlibavaDigitizer::~AllPixAlibavaDigitizer()
{

	rootoutputfile->cd();

	TF1 * fitleft = new TF1("fitleft","landau");
	histochargeleft->Fit(fitleft,"Q");
	histochargeleft->Write();

	TF1 * fitmiddle = new TF1("fitmiddle","landau");
	histochargemiddle->Fit(fitmiddle,"Q");
	histochargemiddle->Write();

	TF1 * fitright = new TF1("fitright","landau");
	histochargeright->Fit(fitright,"Q");
	histochargeright->Write();

	histohitcount->Write();
	histohitcountoverthresh->Write();
	histoclustersize->Write();

	histoloss_position->Write();
	histoloss_share->Write();

	TF1 * fitadc = new TF1("fitadc","landau");
	histoadcout->Fit(fitadc,"Q");
	histoadcout->Write();
	TF1 * fitelec = new TF1("fitelec","landau");
	histoeleout->Fit(fitelec,"Q");
	histoeleout->Write();

	rootoutputfile->Close();

}

void AllPixAlibavaDigitizer::Digitize()
{

	// create the digits collection
	m_digitsCollection = new AllPixAlibavaDigitsCollection("AllPixAlibavaDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;

	// the pixel pointed to
	pair<G4int, G4int> tempPixel;

	// the left neighbour
	pair<G4int, G4int> leftcsharePixel;

	// the right neighbour
	pair<G4int, G4int> rightcsharePixel;

	// the entry count
	G4int nEntries = hitsCollection->entries();

	// summed signals in the 3 channels
	G4double totalcentersignal = 0.0;
	G4double totalleftsignal = 0.0;
	G4double totalrightsignal = 0.0;

	// summed count of the charge lost due to position
	G4double totallosssignal = 0.0;

	// summed count of the charge lost due to charge sharing
	G4double totalshare = 0.0;

	// the output signal in electrons
	G4double electronsignal = 0.0;

	// the output signal in adcs
	G4double adcsignal = 0.0;

	// electrons created per eV
	G4double elec = 3.64;

	// the loop
	for(G4int itr  = 0 ; itr < nEntries ; itr++)
	{

		// positions
		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel()/um;

		// only y for now...

		// the signal for the individual pixels
		G4double leftsignal = 0.0;
		G4double centersignal = 0.0;
		G4double rightsignal = 0.0;

		// the signal in electrons
		G4int integerelectron = (((*hitsCollection)[itr]->GetEdep())/eV)/elec;
		electronsignal += integerelectron;

		// the signal with gain (in MeV)
		//G4double signal = (gain*((*hitsCollection)[itr]->GetEdep()));
		G4double signal = gain * integerelectron;
		adcsignal += signal;

		//G4cout << " " << G4endl;
		//G4cout << "Signal with gain is " << signal << " ADCs" <<G4endl;
		//G4cout << "Actual deposition " << ((*hitsCollection)[itr]->GetEdep())/MeV << " MeV = " << (((*hitsCollection)[itr]->GetEdep())/eV)/elec << " e-" << G4endl;

		// where are we hitting the pixel?
		//G4cout << "Position is " << vec[1] << " um" <<G4endl;

		// reduce charge depending on position relative to the pixel centre
		double positionloss = distancereduce * fabs(vec[1]) * signal;
		//G4cout << "Position loss is " << positionloss << " ADCs" << G4endl;

		// the charge to deposit on the center and to be shared to both neighbours
		double chargetoshare = signal - positionloss;
		//G4cout << "Charge for deposition and sharing is " << chargetoshare << " ADCs" <<G4endl;

		// how much is 'lost' due to the position in the pixel
		totallosssignal += positionloss;
		
		//G4cout << "Charge 'lost' due to position is " << positionloss << " keV" << G4endl;

		// chargeshare spreads a percentage of the actual charge to the two neighbouring pixels
		// this is randomized
		double chargeshare = CLHEP::RandGauss::shoot(sharemean, sharesig);
		//G4cout << "Sharing " << chargeshare*100.0 << " % of the deposition charge" << G4endl;

		// for the center
		centersignal = chargetoshare - (chargetoshare * chargeshare);
		totalshare += chargetoshare * chargeshare;

		// assume this is right...
		if (vec[1] > 0.0)
		{
			rightsignal = (chargetoshare * chargeshare)/2.0;
		}
		if (vec[1] < 0.0)
		{
			leftsignal = (chargetoshare * chargeshare)/2.0;
		}

		// positions left and right
		pixelsContent[tempPixel] +=centersignal;
		leftcsharePixel.first = tempPixel.first;
		leftcsharePixel.second = tempPixel.second - 1;
		rightcsharePixel.first = tempPixel.first;
		rightcsharePixel.second = tempPixel.second + 1;
		pixelsContent[leftcsharePixel] += leftsignal;
		pixelsContent[rightcsharePixel] += rightsignal;

		totalleftsignal += leftsignal;
		totalcentersignal += centersignal;
		totalrightsignal += rightsignal;

		//G4cout << "Left   channel is " << leftcsharePixel.second << " with charge " << pixelsContent[leftcsharePixel] << G4endl;
		//G4cout << "Middle channel is " << tempPixel.second << " with charge " << pixelsContent[tempPixel] << G4endl;
		//G4cout << "Right  channel is " << rightcsharePixel.second << " with charge " << pixelsContent[rightcsharePixel] << G4endl;

	}

	//G4cout << " " << G4endl;
	//G4cout << "Total Track output signal is " << adcsignal << " ADCs, or " << electronsignal << " e-" << G4endl;

	if (adcsignal>0.0)
	{
		histoadcout->Fill(adcsignal);
	}
	if (electronsignal>0.0)
	{
		histoeleout->Fill(electronsignal);
	}
	if (totalleftsignal > 0.0)
	{
		histochargeleft->Fill(totalleftsignal);
	}
	if (totalcentersignal > 0.0)
	{
		histochargemiddle->Fill(totalcentersignal);
	}
	if (totalrightsignal > 0.0)
	{
		histochargeright->Fill(totalrightsignal);
	}
	histoloss_position->Fill(totallosssignal);
	histoloss_share->Fill(totalshare);
	double noise = CLHEP::RandGauss::shoot(noisemean, noisesigma);

	int clustersize = 0;
	if ((totalcentersignal+noise) > seedcut)
	{
		clustersize++;

		noise = CLHEP::RandGauss::shoot(noisemean, noisesigma);
		if ((totalleftsignal+noise) > clustercut)
		{
			clustersize++;
		}
		noise = CLHEP::RandGauss::shoot(noisemean, noisesigma);
		if ((totalrightsignal+noise) > clustercut)
		{
			clustersize++;
		}
	}
	histoclustersize->Fill(clustersize);

	// Now create digits, one per pixel
	// Second entry in the map is the energy deposit in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	// NOTE that there is a nice interface which provides useful info for hits.
	// For instance, the following member gives you the position of a hit with respect
	//  to the center of the pixel.
	// G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
	// See class AllPixTrackerHit !

	// Also, you have full access to the detector geometry from this scope
	// AllPixGeoDsc * GetDetectorGeoDscPtr()
	// provides you with a pointer to the geometry description.
	// See class AllPixGeoDsc !
	G4int hitcount = 0;
	G4int hitcountoverthresh = 0;

	for( ; pCItr != pixelsContent.end() ; pCItr++)
	{
		hitcount++;
		//if((*pCItr).second > m_digitIn.thl) // over threshold !
		if((*pCItr).second > threshold) // over threshold !
		{
			hitcountoverthresh++;
			G4double signal = 0.0;
			AllPixAlibavaDigit * digit = new AllPixAlibavaDigit;
			digit->SetPixelIDX((*pCItr).first.first);
			digit->SetPixelIDY((*pCItr).first.second);
			signal = (*pCItr).second;
			digit->SetPixelCounts(signal);
			m_digitsCollection->insert(digit);
		}
	}

	histohitcount->Fill(hitcount);
	histohitcountoverthresh->Fill(hitcountoverthresh);

	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0)
	{
		G4cout << "--------> Digits Collection : " << collectionName[0] << "(" << m_hitsColName[0] << ")" << " contains " << dc_entries << " digits" << G4endl;
	}

	StoreDigiCollection(m_digitsCollection);

}
