/**
 * Author:
 *    thomas eichhorn <thomas.eichhorn@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixtestalibavaDigitizer_h
#define AllPixtestalibavaDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixtestalibavaDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

using namespace std;

/**
 *  Digitizer AllPixtestalibava implementation
 */
class AllPixtestalibavaDigitizer : public  AllPixDigitizerInterface {

public:

  AllPixtestalibavaDigitizer(G4String, G4String, G4String);
  virtual ~AllPixtestalibavaDigitizer();

  void SetPrimaryVertex(G4PrimaryVertex * pv) {m_primaryVertex = pv;};
  void Digitize ();
  void SetDetectorDigitInputs(G4double){};

private:

  // digitInput typedef is defined in AllPixDigitizerInterface.hh
  digitInput m_digitIn;

  void InitVariables();

	G4double noisemean;
	G4double noisesig;
	G4double pedestalmean;
	G4double pedestalsig;
	G4double commonmodemean;
	G4double commonmodesig;
	G4double sharemean;
	G4double sharesig;
	G4double distancereduce;

	int NPixelX;
	int NPixelY;
	G4double pitchX;
	G4double pitchY;

  AllPixtestalibavaDigitsCollection * m_digitsCollection;
  vector<G4String> m_hitsColName;
  G4PrimaryVertex * m_primaryVertex; // information from EventAction

};

#endif
