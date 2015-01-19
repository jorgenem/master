// main24.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to run SUSY processes in Pythia8.
// All input is specified in the main22.cmnd file.



#include "Pythia8/Pythia.h"
#include <iostream>
#include <vector>

using namespace Pythia8;
using namespace std;




int Sign(double x) {
  if (x >= 0) {
    return 1;
  }
  else {
    return -1;
  }
}
double Invmass(Vec4 p) {
  double invmasssquared =  pow(p.e(),2) - pow(p.px(),2) - pow(p.py(),2) - pow(p.pz(),2);
  return Sign(invmasssquared)*sqrt(abs( invmasssquared ));
}
bool isInAcceptedList (Particle p) {
  bool isInAcceptedList = false;
  int acceptedList[] = {1000001,1000002,1000003,1000004,1000021};
  for (int i = 0; i<5; i++)
  {
    // cout << acceptedList[i] << endl;
    if (abs(p.id()) == acceptedList[i])
    {
      isInAcceptedList = true;
    }
  }
  return isInAcceptedList;
}




int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main24.cmnd");

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors");
  // double eCM   = pythia.parm("Beams:eCM");

  // Initialize.
  pythia.init();

  // Histograms.
  // double epTol = 1e-6 * eCM;
  // Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  // Hist nFinal("final particle multiplicity",100,-0.5,799.5);
  // Hist dnparticledy("dn/dy for particles",100,-10.,10.);

  // Begin event loop.
  int iAbort = 0;
  int iEvent = 0;
  int nCorrect = 0;
  while (iEvent < nEvent) {

    // Generate events. Quit if failure.
    GENERATE: if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    cout << "Event number " << (iEvent+1) << endl;
    ++iEvent;


    // Check that primary hard process is left-handed first- or second generation squarks or gluinos
    bool correctPrimaryProcess = ( isInAcceptedList(event[5]) && isInAcceptedList(event[6]) );
    if (!correctPrimaryProcess) 
    {
      cout << "Break: Incorrect primary process." << endl;
      goto GENERATE;
    }

    bool chain[2];

    for (int iC = 0; iC<2; iC++) // Loop over the two primary outgoing
    {
      if (isInAcceptedList(event[iC+5])) 
      {
        int primaryOutgoing = 5+iC;
        int squark;
        if (abs(event[primaryOutgoing].id()) >= 1000001 && abs(event[primaryOutgoing].id()) <= 1000004)
        {
          // it's a left-handed squark of first or seond gen
          squark = event[primaryOutgoing].iBotCopyId();
        }
        else if (abs(event[primaryOutgoing].id() == 1000021)) 
        { 
          // it's a gluino
          int child1 = event.daughterList(event[5+iC].iBotCopyId())[0];
          int child2 = event.daughterList(event[5+iC].iBotCopyId())[1];
          if (abs(event[child1].id()) >= 1000001 && abs(event[child1].id()) <= 1000004)
          {
            squark = child1;
          }
          else if (abs(event[child2].id()) >= 1000001 && abs(event[child2].id()) <= 1000004)
          {
            squark = child2;
          }
          else 
          {
            cout << "Break: Chain number " << iC+1 << ", was gluino, but no squark children found." << endl;
            cout << "Gluino children have IDs " << event[child1].id() << " and " << event[child2].id() << endl;
            goto GENERATE;
          }
        }

        // We have a squark.
        // Now to check its decay products

        int neutralino2, quark;
        if ( abs(event[event.daughterList(event[squark].iBotCopyId())[0]].id()) == 1000023 
            && 
            ( abs(event[event.daughterList(event[squark].iBotCopyId())[1]].id()) == 1 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[1]].id()) == 2 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[1]].id()) == 3 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[1]].id()) == 4 
            )
          ) 
        {
          neutralino2 = event.daughterList(event[squark].iBotCopyId())[0];
          quark       = event.daughterList(event[squark].iBotCopyId())[1];
        }
        else if ( abs(event[event.daughterList(event[squark].iBotCopyId())[1]].id()) == 1000023 
            && 
            ( abs(event[event.daughterList(event[squark].iBotCopyId())[0]].id()) == 1 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[0]].id()) == 2 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[0]].id()) == 3 || 
              abs(event[event.daughterList(event[squark].iBotCopyId())[0]].id()) == 4 
            ) 
          )
        {
          neutralino2 = event.daughterList(event[squark].iBotCopyId())[1];
          quark       = event.daughterList(event[squark].iBotCopyId())[0];
        }
        else
        {
          cout << "Break: Chain number " << iC+1 << ", squark didn't have the right children." << endl;
          goto GENERATE;
        }

        // We have a neutralino2 and a quark
        // Now to check neutralino2's decay products

        int slepton, lepton1;
        if  (  
              (
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 1000011 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 1000013 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 2000011 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 2000013
              )
              &&
              (
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 11 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 13
              )
            )
        {
          slepton = event.daughterList(event[neutralino2].iBotCopyId())[0];
          lepton1 = event.daughterList(event[neutralino2].iBotCopyId())[1];
        }
        else if (
              (
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 1000011 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 1000013 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 2000011 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[1]].id()) == 2000013
              )
              &&
              (
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 11 ||
                abs(event[event.daughterList(event[neutralino2].iBotCopyId())[0]].id()) == 13
              )
             )
        {
          slepton = event.daughterList(event[neutralino2].iBotCopyId())[1];
          lepton1 = event.daughterList(event[neutralino2].iBotCopyId())[0];
        }
        else
        {
          cout << "Break: Chain number " << iC+1 << ", neutralino2 didn't have the right children." << endl;
          goto GENERATE;
        }

        // We have a slepton and a lepton
        // Now to check slepton's decay products

        int neutralino1, lepton2;
        if  (  
              (
                abs(event[event.daughterList(event[slepton].iBotCopyId())[0]].id()) == 1000022
              )
              &&
              (
                abs(event[event.daughterList(event[slepton].iBotCopyId())[1]].id()) == 11 ||
                abs(event[event.daughterList(event[slepton].iBotCopyId())[1]].id()) == 13
              )
            )
        {
          neutralino1 = event.daughterList(event[slepton].iBotCopyId())[0];
          lepton2     = event.daughterList(event[slepton].iBotCopyId())[1];
        }
        else if (
              (
                abs(event[event.daughterList(event[slepton].iBotCopyId())[1]].id()) == 1000022
              )
              &&
              (
                abs(event[event.daughterList(event[slepton].iBotCopyId())[0]].id()) == 11 ||
                abs(event[event.daughterList(event[slepton].iBotCopyId())[0]].id()) == 13
              )
             )
        {
          neutralino1 = event.daughterList(event[slepton].iBotCopyId())[1];
          lepton2     = event.daughterList(event[slepton].iBotCopyId())[0];
        }
        else
        {
          cout << "Break: Chain number " << iC+1 << ", slepton didn't have the right children." << endl;
          goto GENERATE;
        }

        // Now we have all the correct particles in the chain.
        chain[iC] = true;

        // Print them to screen for testing
        cout << "Correct chain instance found." << endl;
        cout << "Squark  ID = " << event[squark].id() << endl;
        cout << "Chi2    ID = " << event[neutralino2].id() << endl;
        cout << "Quark   ID = " << event[quark].id() << endl;
        cout << "Slepton ID = " << event[slepton].id() << endl;
        cout << "Lepton1 ID = " << event[lepton1].id() << endl;
        cout << "Chi1    ID = " << event[neutralino1].id() << endl;
        cout << "Lepton2 ID = " << event[lepton2].id() << endl;

        } // END IF it is a squark/gluino
        else 
        {
          cout << "Break: Chain number " << iC+1 << " didn't have a squark or gluino" << endl;
          goto GENERATE;
        }


      } // END FOR loop over chains

      // If we have made it this far, both chains are correct
      ++nCorrect;
      

      

    //   // cout << event[5+iC].id() << " " << event[event.daughterList(event[5+iC].iBotCopy())[0]].id() << " " << event[event.daughterList(event[5+iC].iBotCopy())[1]].id() << endl;   
    // }
    // if (correctPrimaryProcess)
    // {
    //   cout << "isInAcceptedList = true" << ", id = " << event[5].id() << " and " << event[6].id() << endl;
    // }
    // else
    // {
    //   cout << "isInAcceptedList = false" << ", id = " << event[5].id() << " and " << event[6].id() << endl;
    // }



    // // Investigate the primary hard process
    // cout << "Information about primary hard process: \t" << endl;
    // cout << "ID's of primary hard process participants: " << pythia.info.id1() << " \t " << pythia.info.id2() << endl;
    // cout << "Size of primary hard process: " << pythia.process.size() << endl;
    // cout << "ID of event listing 3 and 4: " << event[3].id() << ",\t " << event[4].id() << endl;
    // cout << "Daughterlist of 3: " << event.daughterList(3)[0] << " \t" << event.daughterList(3)[1] << endl;
    // cout << "ID of event listing 5 and 6: " << event[5].id() << ",\t " << event[6].id() << endl;
    // cout << "Daughterlist of 5: " << event.daughterList(5)[0] << " \t" << event.daughterList(6)[1] << endl;
    // cout << "ID, daughterlist of 9: " << event[9].id() << " \t" << event.daughterList(9)[0] << " \t" << event.daughterList(9)[1] << endl;
    // cout << "Bottom carbon copy of 5: " << event[5].iBotCopy() /*<< " \t" << event[event.iBotCopy(5)].id()*/ << endl;
    // cout << "4-momentum and invariant mass of 9: \t" << event[9].px() << " \t" << event[9].py() << " \t" << event[9].pz() << " \t" << event[9].e() << " \t" << event[9].p().mCalc() << " \t" << Invmass(event[9].p()) << endl;
    // cout << "Daughterlist and ID's of 9 : \t" << event.daughterList(9)[0] << " \t" << event[event.daughterList(9)[0]].id() << " \t" << event.daughterList(9)[1] << " \t" << event[event.daughterList(9)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 12: \t" << event[12].px() << " \t" << event[12].py() << " \t" << event[12].pz() << " \t" << event[12].e() << " \t" << event[12].p().mCalc() << " \t" << Invmass(event[12].p()) << endl;
    // cout << "Daughterlist and ID's of 12 : \t" << event.daughterList(12)[0] << " \t" << event[event.daughterList(12)[0]].id() << " \t" << event.daughterList(12)[1] << " \t" << event[event.daughterList(12)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 20: \t" << event[20].px() << " \t" << event[20].py() << " \t" << event[20].pz() << " \t" << event[20].e() << " \t" << event[20].p().mCalc() << " \t" << Invmass(event[20].p()) << endl;
    // cout << "Daughterlist and ID's of 20 : \t" << event.daughterList(20)[0] << " \t" << event[event.daughterList(20)[0]].id() << " \t" << event.daughterList(20)[1] << " \t" << event[event.daughterList(20)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 38: \t" << event[38].px() << " \t" << event[38].py() << " \t" << event[38].pz() << " \t" << event[38].e() << " \t" << event[38].p().mCalc() << " \t" << Invmass(event[38].p()) << endl;
    // cout << "Daughterlist and ID's of 38 : \t" << event.daughterList(38)[0] << " \t" << event[event.daughterList(38)[0]].id() << " \t" << event.daughterList(38)[1] << " \t" << event[event.daughterList(38)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 74: \t" << event[74].px() << " \t" << event[74].py() << " \t" << event[74].pz() << " \t" << event[74].e() << " \t" << event[74].p().mCalc() << " \t" << Invmass(event[74].p()) << endl;
    // cout << "Daughterlist and ID's of 74 : \t" << event.daughterList(74)[0] << " \t" << event[event.daughterList(74)[0]].id() << " \t" << event.daughterList(74)[1] << " \t" << event[event.daughterList(74)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 152: \t" << event[152].px() << " \t" << event[152].py() << " \t" << event[152].pz() << " \t" << event[152].e() << " \t" << event[152].p().mCalc() << " \t" << Invmass(event[152].p()) << endl;
    // cout << "Daughterlist and ID's of 152 : \t" << event.daughterList(152)[0] << " \t" << event[event.daughterList(152)[0]].id() << " \t" << event.daughterList(152)[1] << " \t" << event[event.daughterList(152)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 199: \t" << event[199].px() << " \t" << event[199].py() << " \t" << event[199].pz() << " \t" << event[199].e() << " \t" << event[199].p().mCalc() << " \t" << Invmass(event[199].p()) << endl;
    // cout << "Daughterlist and ID's of 199 : \t" << event.daughterList(199)[0] << " \t" << event[event.daughterList(199)[0]].id() << " \t" << event.daughterList(199)[1] << " \t" << event[event.daughterList(199)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 200: \t" << event[200].px() << " \t" << event[200].py() << " \t" << event[200].pz() << " \t" << event[200].e() << " \t" << event[200].p().mCalc() << " \t" << Invmass(event[200].p()) << endl;
    // cout << "Daughterlist and ID's of 200 : \t" << event.daughterList(200)[0] << " \t" << event[event.daughterList(200)[0]].id() << " \t" << event.daughterList(200)[1] << " \t" << event[event.daughterList(200)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 201: \t" << event[201].px() << " \t" << event[201].py() << " \t" << event[201].pz() << " \t" << event[201].e() << " \t" << event[201].p().mCalc() << " \t" << Invmass(event[201].p()) << endl;
    // cout << "Daughterlist and ID's of 201 : \t" << event.daughterList(201)[0] << " \t" << event[event.daughterList(201)[0]].id() << " \t" << event.daughterList(201)[1] << " \t" << event[event.daughterList(201)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 206: \t" << event[206].px() << " \t" << event[206].py() << " \t" << event[206].pz() << " \t" << event[206].e() << " \t" << event[206].p().mCalc() << " \t" << Invmass(event[206].p()) << endl;
    // cout << "Daughterlist and ID's of 206 : \t" << event.daughterList(206)[0] << " \t" << event[event.daughterList(206)[0]].id() << " \t" << event.daughterList(206)[1] << " \t" << event[event.daughterList(206)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 207: \t" << event[207].px() << " \t" << event[207].py() << " \t" << event[207].pz() << " \t" << event[207].e() << " \t" << event[207].p().mCalc() << " \t" << Invmass(event[207].p()) << endl;
    // cout << "Daughterlist and ID's of 207 : \t" << event.daughterList(207)[0] << " \t" << event[event.daughterList(207)[0]].id() << " \t" << event.daughterList(207)[1] << " \t" << event[event.daughterList(207)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 212: \t" << event[212].px() << " \t" << event[212].py() << " \t" << event[212].pz() << " \t" << event[212].e() << " \t" << event[212].p().mCalc() << " \t" << Invmass(event[212].p()) << endl;
    // cout << "Daughterlist and ID's of 212 : \t" << event.daughterList(212)[0] << " \t" << event[event.daughterList(212)[0]].id() << " \t" << event.daughterList(212)[1] << " \t" << event[event.daughterList(212)[1]].id() << endl; 
    // cout << "4-momentum and invariant mass of 215: \t" << event[215].px() << " \t" << event[215].py() << " \t" << event[215].pz() << " \t" << event[215].e() << " \t" << event[215].p().mCalc() << " \t" << Invmass(event[215].p()) << endl;
    // cout << "Daughterlist and ID's of 215 : \t" << event.daughterList(215)[0] << " \t" << event[event.daughterList(215)[0]].id() << " \t" << event.daughterList(215)[1] << " \t" << event[event.daughterList(215)[1]].id() << endl; 

    // Loop over all (prev:final) particles in the event.
    // cout << "Information about all hadrons in event: " << endl;
    int nFin = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) /*if (event[i].isFinal())*/ {
      nFin++;
      pSum += event[i].p();
      // dnparticledy.fill(event[i].y());
        if (abs(event[i].id())<10) {
        // cout << "ID \t INVMASS \t LISTMASS" << endl;
        // double invmasssquared = pow(event[i].e(),2) - pow(event[i].px(),2) - pow(event[i].py(),2) - pow(event[i].pz(),2);
        // cout << event[i].id() << " \t " << Sign(invmasssquared) << "*" << sqrt(abs(invmasssquared)) << " \t " << event[i].m() << endl;
      }
    }

    // Check and print event with too big energy-momentum deviation.
    // nFinal.fill(nFin);
    // double epDev = abs(pSum.e() - eCM) + abs(pSum.px()) + abs(pSum.py())
    //   + abs(pSum.pz());
    // epCons.fill(epDev);
    // if (epDev > epTol) {
    //   cout << " Warning! Event with epDev = " << scientific
    //        << setprecision(4) << epDev << " now listed:";
    //   event.list();
    // }

  // End of event loop.
  }

  // Final statistics and histogram output.
  // pythia.stat();
  // cout << epCons << nFinal << dnparticledy;

  return 0;
}

