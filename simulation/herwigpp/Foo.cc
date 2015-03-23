// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Foo class.
//

#include "Foo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <armadillo>
#include <cmath>



using namespace Jorgen;
//using namespace arma;

Foo::Foo() {}

Foo::~Foo() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

// ParticleVector getChildren(PPtr p)
// {
//   int id = p->id();
//   ParticleVector children = p->children();
//   if(children.size() > 0 && children[0]->id() == id)
//   {
//     return getChildren(children[0]);
//   }
//   else if(children.size() == 0)
//   {
//     cout << "Doesn't have children!" << endl;
//   }
//   else
//   {
//     return children;
//   }
// }

bool getChildren2(PPtr p, ParticleVector &children)
{
  int id = p->id();
  children = p->children();
  bool hasChildren = (children.size() > 0);
  if(!hasChildren)
  {
    return false;
  }
  if (children[0]->id() == id)
  {
    return getChildren2(children[0], children);
  }
  else
  {
    return true;
  }
}

PPtr getFinalInstance(PPtr p) // This function works its way down in the instance chain until p's first child is not p, or p has no children.
{
  //int id = p->id();
  ParticleVector children = p->children();
  if(children.size()==0) 
  {
    return p;
  }
  else if (children[0]->id() != p->id()) // Is this a possible flaw? Are we certain that the particle itself is always indexed first if it e.g. emits gluons or gammas?
  {
    return p;
  }
  else
  {
    return getFinalInstance(children[0]);
  }
}
PPtr getFinalInstanceBeforeRadiation(PPtr p) // Added 20150116 as a variation: The point is to return p if p has more than one child and first child is itself (i.e. if it radiates something in the next step)
{
  //int id = p->id();
  ParticleVector children = p->children();
  if(children.size()==0) 
  {
    return p;
  }
  else if (children[0]->id() == p->id() && children.size()>1) // added 20150116
  {
    return p;
  }
  else if (children[0]->id() != p->id()) // Is this a possible flaw? Are we certain that the particle itself is always indexed first if it e.g. emits gluons or gammas?
  {
    return p;
  }
  else
  {
    return getFinalInstance(children[0]);
  }
}

PPtr getOffshellQuark(PPtr p)
{
  ParticleVector children = p->children();
  if (children.size()>1 && children[0]->id() == p->id())
  {
    return p;
  }
  else if (children.size()>1 && children[1]->id() == p->id())
  {
    return p;
  }
  else if (children.size()>1 && children[0]->id() != p->id() && children[1]->id() != p->id() && abs(p->id())<10)
  {
    cout << "WARNING! A quark seems to have decayed. Quark id = " << p->id() << ", children id = " << children[0]->id() << children[1]->id() << ". Check carefully." << endl;
    return p;
  }
  else if (children.size()>1 && children[0]->id() != p->id() && children[1]->id() != p->id())
  {
      // This case applies e.g. for squarks, which should decay.
      return p;
  }
  else if (children.size()==0)
  {
    // cout << "final state quark! id = " << p->id() << endl;
    return p;
  }
  else
  {
    return getOffshellQuark(children[0]);
  }
}

bool momentumConservation(PPtr p1, PPtr p2, PPtr p3)
{
  // Particle p1 decays to particle p2 & p3
  // bool momentumIsConserved = true;
  double epsilon = 0.5; // Should probably not be smaller than BasicConsistency:: maximum 4-momentum violation.
  double p1x = p1->momentum().x();
  double p1y = p1->momentum().y();
  double p1z = p1->momentum().z();
  double p1e = p1->momentum().e();
  double p2x = p2->momentum().x();
  double p2y = p2->momentum().y();
  double p2z = p2->momentum().z();
  double p2e = p2->momentum().e();
  double p3x = p3->momentum().x();
  double p3y = p3->momentum().y();
  double p3z = p3->momentum().z();
  double p3e = p3->momentum().e();

  if (abs(p1x - p2x - p3x)<epsilon
     && abs(p1y - p2y - p3y)<epsilon 
     && abs(p1z - p2z - p3z)<epsilon 
     && abs(p1e - p2e - p3e)<epsilon)
  {
    return true;
  }
  else
  {
    cout << "Momentum conservation failed without gammas: Dpx = " << p1x - p2x - p3x << ", Dpy = " << p1y - p2y - p3y << ", Dpz = " << p1z - p2z - p3z << ", De = " << p1e - p2e - p3e << endl;
    cout << "Individual particle momenta: " << endl;
    cout << "Incoming: \t id = " << p1->id() << " \t\t" << p1x << " \t\t" << p1y << " \t\t" << p1z << " \t\t" << p1e << endl;
      cout << "Outgoing 1: \t id = " << p2->id() << " \t\t" << p2x << " \t\t" << p2y << " \t\t" << p2z << " \t\t" << p2e << endl;
      cout << "Outgoing 2: \t id = " << p3->id() << " \t\t" << p3x << " \t\t" << p3y << " \t\t" << p3z << " \t\t" << p3e << endl;


    return false;
  }
}

bool momentumConservation1to4(PPtr p1, PPtr p2, PPtr p3, PPtr p4, PPtr p5)
{
  // Made 20150314
  // Particle p1 decays to particle p2, p3, p4 and p5
  // bool momentumIsConserved = true;
  double epsilon = 10.0; // NB! MeV unit. Should probably not be smaller than BasicConsistency:: maximum 4-momentum violation.
  double p1x = p1->momentum().x();
  double p1y = p1->momentum().y();
  double p1z = p1->momentum().z();
  double p1e = p1->momentum().e();
  double p2x = p2->momentum().x();
  double p2y = p2->momentum().y();
  double p2z = p2->momentum().z();
  double p2e = p2->momentum().e();
  double p3x = p3->momentum().x();
  double p3y = p3->momentum().y();
  double p3z = p3->momentum().z();
  double p3e = p3->momentum().e();
  double p4x = p4->momentum().x();
  double p4y = p4->momentum().y();
  double p4z = p4->momentum().z();
  double p4e = p4->momentum().e();
  double p5x = p5->momentum().x();
  double p5y = p5->momentum().y();
  double p5z = p5->momentum().z();
  double p5e = p5->momentum().e();

  if (abs(p1x - p2x - p3x - p4x -p5x)<epsilon
     && abs(p1y - p2y - p3y - p4y - p5y)<epsilon 
     && abs(p1z - p2z - p3z - p4z - p5z)<epsilon 
     && abs(p1e - p2e - p3e - p4e - p5e)<epsilon)
  {
    return true;
  }
  else
  {
    // cout << "Momentum conservation failed without gammas: Dpx = " << p1x - p2x - p3x << ", Dpy = " << p1y - p2y - p3y << ", Dpz = " << p1z - p2z - p3z << ", De = " << p1e - p2e - p3e << endl;
    // cout << "Individual particle momenta: " << endl;
    // cout << "Incoming: \t id = " << p1->id() << " \t\t" << p1x << " \t\t" << p1y << " \t\t" << p1z << " \t\t" << p1e << endl;
    //   cout << "Outgoing 1: \t id = " << p2->id() << " \t\t" << p2x << " \t\t" << p2y << " \t\t" << p2z << " \t\t" << p2e << endl;
    //   cout << "Outgoing 2: \t id = " << p3->id() << " \t\t" << p3x << " \t\t" << p3y << " \t\t" << p3z << " \t\t" << p3e << endl;


    return false;
  }
}

bool momentumConservationWithPossibleGammas(PPtr p1, PPtr p2, PPtr p3) {
  // if (abs(p1->id())>=1000001 && abs(p1->id())<=1000006) // If it is a squark, make sure it has correct off-shell momentum
  // {
  //   p1 = getOffshellQuark(p1);
  // } // We should NOT pick the off-shell version of the squark -- we need to pick it out AFTER it has radiated gluons, just before it actually decays to neutralino2+quark. Then we need to make sure we pick the neutralino2 and quark BEFORE they radiate anything.
  if (p1->children().size()==2) 
  {
    return momentumConservation(p1, p2, p3);
  }
  else if (p1->children().size()>2) 
  {
    // cout << "An event with gammas! There are " << p1->children().size()-2 << " of them, the first one has ID = " << p1->children()[2]->id() << endl;
    double epsilon = 5e-1;
    double p1x = p1->momentum().x();
    double p1y = p1->momentum().y();
    double p1z = p1->momentum().z();
    double p1e = p1->momentum().e();
    double p2x = p2->momentum().x();
    double p2y = p2->momentum().y();
    double p2z = p2->momentum().z();
    double p2e = p2->momentum().e();
    double p3x = p3->momentum().x();
    double p3y = p3->momentum().y();
    double p3z = p3->momentum().z();
    double p3e = p3->momentum().e();

    double outgoingSum_x = p2x + p3x;
    double outgoingSum_y = p2y + p3y;
    double outgoingSum_z = p2z + p3z;
    double outgoingSum_e = p2e + p3e;

    for (ParticleVector::const_iterator gamma = p1->children().begin() + 2; gamma != p1->children().end(); ++gamma ) 
    {
      if ((*gamma)->id() != 22) {
        cout << "Found a child which is not a photon! It has id = " << (*gamma)->id() << endl;
        return false;
        break;
      }
      outgoingSum_x = outgoingSum_x + getFinalInstance(*gamma)->momentum().x();
      outgoingSum_y = outgoingSum_y + getFinalInstance(*gamma)->momentum().y();
      outgoingSum_z = outgoingSum_z + getFinalInstance(*gamma)->momentum().z();
      outgoingSum_e = outgoingSum_e + getFinalInstance(*gamma)->momentum().e();
    }
    if (  abs(p1x - outgoingSum_x)<epsilon 
        &&abs(p1y - outgoingSum_y)<epsilon 
        &&abs(p1z - outgoingSum_z)<epsilon 
        &&abs(p1e - outgoingSum_e)<epsilon ) {
      return true;
    }
    else {
      cout << "Momentum conservation failed with gammas included: Dpx = " << p1x - outgoingSum_x << ", Dpy = " << p1y - outgoingSum_y << ", Dpz = " << p1z - outgoingSum_z << ", De = " << p1e - outgoingSum_e << endl;
      cout << "Individual particle momenta: " << endl;
      cout << "Incoming: \t id = " << p1->id() << " \t\t" << p1x << " \t\t" << p1y << " \t\t" << p1z << " \t\t" << p1e << endl;
      int counter = 1;
      for (ParticleVector::const_iterator gamma = p1->children().begin(); gamma != p1->children().end(); ++gamma) {
        cout << "Outgoing " << counter << ": \t id = " << (*gamma)->id() << " \t\t" << getFinalInstance(*gamma)->momentum().x() << " \t\t" << getFinalInstance(*gamma)->momentum().y() << " \t\t" << getFinalInstance(*gamma)->momentum().z() << " \t\t" << getFinalInstance(*gamma)->momentum().e() << endl;
          ++counter;
      }
      return false;
    } 
  } // end else if there are at least three particles. If none of the above is true, something is wrong (less than two children!)
  else {
    cout << "Error in momentumConservationWithPossibleGammas: p1 has less than two children?" << endl;
    if (p1->children().size()==1) {
      cout << "p1 id = " << p1->id() << ", child 1 id = " << p1->children()[0]->id() << endl;
    }
    else {
      cout << "No children. Final state. p1 id = " << p1->id() << endl;
    }
    return false;
  }
}

bool getChildrenWithPossibleGammas (PPtr p, int id11, int id12, int id13, int id14, int id21, int id22, PPtr &child1, PPtr &child2)
{
  // Possible TODO in this function: Return false if magnitude of gamma energy is too large

  int indexId1, indexId2;
  bool child1IsPresent = false;
  bool child2IsPresent = false;
  ParticleVector children;
  if (getChildren2(p, children))
  {

  }
  else
  {
    // If getChildren2 returns false then we should also return false, it means particle has no children
    return false;
  }
  int counter = 0;
  for (ParticleVector::const_iterator child = children.begin(); child != children.end(); ++child)
  {
    // Check that there are no other children than we expect
    if (abs((*child)->id()) == id11 || abs((*child)->id()) == id12 || abs((*child)->id()) == id13 || abs((*child)->id()) == id14)
    {
      // If child is equal to child1 (slepton) then
      child1IsPresent = true;
      indexId1 = counter; 
    }
    if (abs((*child)->id()) == id21 || abs((*child)->id()) == id22)
    {
      // If child is equal to child2 (lepton) then
      child2IsPresent = true;
      indexId2 = counter;
    }
    if  (
      abs((*child)->id()) != 22 
      && (abs((*child)->id()) != id11 && abs((*child)->id()) != id12 && abs((*child)->id()) != id13 && abs((*child)->id()) != id14) 
      && (abs((*child)->id()) != id21 && abs((*child)->id()) != id22)
        )
    {
      // If child is none of the particles we expect, then break and return false
      return false;
    }
    ++counter;
  }
  if (child1IsPresent && child2IsPresent)
  {
    child1 = getFinalInstance(children[indexId1]);
    child2 = getFinalInstance(children[indexId2]);
    return true;
  }
  else
  {
    return false;
  }



}

struct MomentumVector
{
  int id;
  double m;
  double e;
  double px;
  double py;
  double pz;
};

void Foo::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).

  if ( loop > 0 || state != 0 || !event ) return;
  numEventsTotal++; // Count event as one of the total regardless of whether it is the right type
  cout << numEventsTotal << endl;

  /** Get the final-state particles */
  tPVector particles=event->getFinalState();

  /** Make and open text output file */
  ofstream textOutput;
  // ofstream adOut;
  textOutput.open("LHC-MSSM-analysis.log", ios::out | ios::app);

  // Declare variables 
  vector<MomentumVector> event_particles; // Save all chain particles
  bool chain = true; // Assume that both chains are present until proven wrong by if-tests
  bool momentumIsConserved = true; // Assume that momentum conservation tests succeed in both chains until proven wrong
  bool gluino = false; // Assume we didn't produce a gluino in the hard subprocess
  // Get primary subprocess (squarks):
  tSubProPtr primarysubprocess = event->primarySubProcess();
  ParticleVector squarks = primarysubprocess->outgoing();

  // Check that we produced two left-handed squarks of first
  // or second generation
  // textOutput << "==== Total event number " << numEventsTotal << " ====" << endl;

  // textOutput << squarks[0]->id() << squarks[1]->id() << endl;
  // textOutput << "squarktest = " << (abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004) << squarks.size() << endl;
  ParticleVector gluinochildren; // Declare this in case one of the outgoing hard partons is a gluino

  if (
    (abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004) // both are squarks
    ||
    ( (abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) == 1000021 && getChildren2(squarks[1], gluinochildren)) || (abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004 && abs(squarks[0]->id()) == 1000021 && getChildren2(squarks[0], gluinochildren)) )
    )

  {


    // If one of the hard partons is a gluino, check that the squark child is left-handed first- or second gen
    if ( (abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) == 1000021 && getChildren2(squarks[1], gluinochildren)) || (abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004 && abs(squarks[0]->id()) == 1000021 && getChildren2(squarks[0], gluinochildren)) )
    {
      if (!( 
            (abs(gluinochildren[0]->id()) >= 1000001 && abs(gluinochildren[0]->id()) <= 1000004)
          ||(abs(gluinochildren[1]->id()) >= 1000001 && abs(gluinochildren[1]->id()) <= 1000004)
          ))
        return;
    }


    // cout << "loop over squarks" << endl;
    // Loop over squarks
    for (ParticleVector::const_iterator squarkiterator = squarks.begin(); squarkiterator != squarks.end(); ++squarkiterator)
    {
      // One of the two squarks might be a gluino. Need to check:
      PPtr squark; // this will be the squark, either directly or as a gluino child.
      // If there is a gluino:
      if ((*squarkiterator)->id() == 1000021 && getChildren2(*squarkiterator, gluinochildren))
      {
        // Check which of the children is the squark
        if (abs(gluinochildren[0]->id()) > abs(gluinochildren[1]->id()))
        {
          squark = gluinochildren[0];
        }
        else
        {
          squark = gluinochildren[1];
        }
      }
      else // else it is a squark
      {
        squark = *squarkiterator;
      }
      // cout << "entered squark loop" << endl;
      // Get squark children -- we know there have to be some, since it's a squark
      ParticleVector squarkchildren;
      if (getChildren2(squark, squarkchildren))
      {
        // cout << "getchildren2(squarks)" << endl;
        if ( squarkchildren.size()==2 && ( (abs(squarkchildren[0]->id()) >=1 && abs(squarkchildren[0]->id()) <=4 && squarkchildren[1]->id() == 1000023) || (abs(squarkchildren[1]->id()) >=1 && abs(squarkchildren[1]->id()) <=4 && squarkchildren[0]->id() == 1000023) ) )
        {
          // cout << "squark has the right children" << endl;
          // Check which child is which, make sure the quark is showered
          PPtr neutralino2, quark;
          if (abs(squarkchildren[0]->id()) > abs(squarkchildren[1]->id()))
          {
          //  cout << "trying to get offshell quark of id " << squarkchildren[1]->id() << endl;
            neutralino2 = squarkchildren[0];
            quark = getOffshellQuark(squarkchildren[1]);
          }
          else
          {
           // cout << "trying to get offshell quark of id " << squarkchildren[0]->id() << endl;
            neutralino2 = squarkchildren[1];
            quark = getOffshellQuark(squarkchildren[0]);
          }

          // cout << "got both squark children" << endl;
          // Get neutralino2 children -- has to have children since we know it's a neutralino2
          // ParticleVector neutralino2children;
          PPtr slepton, lepton1;
          if (getChildrenWithPossibleGammas(neutralino2,1000011,1000013,2000011,2000013,11,13, slepton, lepton1))
          {
            // cout << "Got neutralino2 children with possible gammas" << endl;
            if ( true
              // (neutralino2children.size()==2 || neutralino2DecayedThreeBodyWithGamma)
              // && 
              // ( ( (abs(neutralino2children[0]->id()) ==11 || abs(neutralino2children[0]->id()) ==13) && ( abs(neutralino2children[1]->id()) == 1000011 || abs(neutralino2children[1]->id()) == 1000013 || abs(neutralino2children[1]->id()) == 2000011 || abs(neutralino2children[1]->id()) == 2000013 ) ) || ( (abs(neutralino2children[1]->id()) ==11 || abs(neutralino2children[1]->id()) ==13) && ( abs(neutralino2children[0]->id()) == 1000011 || abs(neutralino2children[0]->id()) == 1000013 || abs(neutralino2children[0]->id()) == 2000011 || abs(neutralino2children[0]->id()) == 2000013 ) ) ) 
              )
            {
              // cout << "neutralino2 has the right children" << endl;
              // Check which child is which
              // if (abs(neutralino2children[0]->id()) > abs(neutralino2children[1]->id()))
              // {
              //   slepton = neutralino2children[0];
              //   lepton1 = neutralino2children[1];
              // }
              // else
              // {
              //   slepton = neutralino2children[1];
              //   lepton1 = neutralino2children[0];      
              // }

              // cout << "got both neutralino2 children" << endl;
              // Get slepton children -- has to have children since we know it's a slepton
              // ParticleVector sleptonchildren;>
              PPtr neutralino1, lepton2;
              if (getChildrenWithPossibleGammas(slepton, 1000022,1000022,1000022,1000022,11,13, neutralino1, lepton2))
              {
                // cout << "Got slepton children with possible gammas" << endl;
                // bool sleptonDecayedThreeBodyWithGamma = (sleptonchildren.size()==3 && sleptonchildren[2]->id()==22);
                if ( true
                  // (sleptonchildren.size()==2 || sleptonDecayedThreeBodyWithGamma)
                  // && 
                  // ( ( (abs(sleptonchildren[0]->id()) ==11 || abs(sleptonchildren[0]->id()) ==13) && abs(sleptonchildren[1]->id()) == 1000022 ) || ( (abs(sleptonchildren[1]->id()) ==11 || abs(sleptonchildren[1]->id()) ==13) && abs(sleptonchildren[0]->id()) == 1000022 ) ) 
                  )
                {
                  // cout << "slepton has the right children" << endl;
                  // PPtr neutralino1, lepton2;
                  // if (abs(sleptonchildren[0]->id()) > abs(sleptonchildren[1]->id()))
                  // {
                  //   neutralino1 = sleptonchildren[0];
                  //   lepton2 = sleptonchildren[1];
                  // }
                  // else
                  // {
                  //   neutralino1 = sleptonchildren[1];
                  //   lepton2 = sleptonchildren[0];
                  // }

                  // cout << "got both slepton children" << endl;



                  // Check momentum conservation, accepting events with radiation of photons that are lost
                  // PPtr squark = getOffshellQuark(*squarkiterator);
                  // PPtr squarkoffshell = getOffshellQuark(squark);
                  PPtr squarkfinal = getFinalInstance(squark);
                  PPtr neutralino2final = getFinalInstance(neutralino2);
                  PPtr neutralino2beforeradiation = getFinalInstanceBeforeRadiation(neutralino2);
                  PPtr lepton1final = getFinalInstance(lepton1);
                  PPtr sleptonfinal = getFinalInstance(slepton);
                  PPtr sleptonbeforeradiation = getFinalInstanceBeforeRadiation(slepton);
                  PPtr lepton2final = getFinalInstance(lepton2);
                  PPtr neutralino1beforeradiation = getFinalInstanceBeforeRadiation(neutralino1);


                  bool squarkDecayMomentumConservation = momentumConservationWithPossibleGammas(squarkfinal, neutralino2beforeradiation, quark);
                  bool neutralino2finalDecayMomentumConservation = momentumConservationWithPossibleGammas(neutralino2final, sleptonbeforeradiation, lepton1final);
                  bool sleptonDecayMomentumConservation = momentumConservationWithPossibleGammas(sleptonfinal, lepton2final, neutralino1beforeradiation);
                  if (!(neutralino2finalDecayMomentumConservation && sleptonDecayMomentumConservation && squarkDecayMomentumConservation)) 
                  {
                    return;
                  }




                  // if (!(neutralino2finalDecayMomentumConservation && sleptonDecayMomentumConservation && squarkDecayMomentumConservation) &&  abs(squark->id())>=1000001 && abs(squark->id())<=1000004) // Check that it's not a chain with right-handed squarks)
                  // {
                  //   // textOutput << "=== Total event number " << numEventsTotal << " (failed momentum conservation) ===" << endl;
                  //   // textOutput << "Momentum conservation test failed." << endl;
                  //   // textOutput << "squark conserved = " << squarkDecayMomentumConservation << ", \t id = " << squarkfinal->id() << endl;
                  //   // textOutput << "squark momentum \t" << squarkfinal->momentum().x() << "\t" << squarkfinal->momentum().y() << "\t" << squarkfinal->momentum().z() << "\t" << squarkfinal->momentum().e() << "\t" <<
                  //   //   squarkfinal->momentum().m() << endl;
                  //   // textOutput << "neutralino2final momentum \t" << neutralino2final->momentum().x() << "\t" << neutralino2final->momentum().y() << "\t" << neutralino2final->momentum().z() << "\t" << neutralino2final->momentum().e() << "\t" <<
                  //   //   neutralino2final->momentum().m() << endl;
                  //   // textOutput << "quark momentum \t" << quark->momentum().x() << "\t" << quark->momentum().y() << "\t" << quark->momentum().z() << "\t" << quark->momentum().e() << "\t" <<
                  //   //   quark->momentum().m() << endl;
                  //   // textOutput << "neutralino2final conserved = " << neutralino2finalDecayMomentumConservation << ", \t id = " << neutralino2final->id() << endl;
                  //   // // textOutput << "neutralino2 decayed three-body with gamma = " << neutralino2DecayedThreeBodyWithGamma << endl;
                  //   // textOutput << "slepton conserved = " << sleptonDecayMomentumConservation << ", \t id = " << sleptonfinal->id() << endl;
                  //   // // textOutput << "slepton decayed three-body with gamma = " << sleptonDecayedThreeBodyWithGamma << endl;

                  //   momentumIsConserved = false;
                  // }
                  // else
                  // {
                  //   // textOutput << "Momentum conservation test succeeded." << endl;
                  //   // textOutput << "squark conserved = " << squarkDecayMomentumConservation << endl;
                  //   // textOutput << "neutralino2final conserved = " << neutralino2finalDecayMomentumConservation << endl;
                  //   // textOutput << "slepton conserved = " << sleptonDecayMomentumConservation << endl;
                  //   // numEventsMomentumConservation++;
                  // }



                  // 20150314 Check momentum conservation of the four end products versus the squark. Disregard all events where it is not conserved,
                  // this probably means any event with gammas
                  // squark = getFinalInstance(squark);
                  // bool momentumIsConserved = momentumConservation1to4(squark, quark, lepton1, lepton2, neutralino1);
                  // if (!momentumIsConserved)
                  //   return;


                  // DEBUG: Save squark 4-momentum
                  MomentumVector Psquark;
                  Psquark.id = squark->id();
                  Psquark.m = squark->momentum().m();
                  Psquark.e = squark->momentum().e();
                  Psquark.px = squark->momentum().x();
                  Psquark.py = squark->momentum().y();
                  Psquark.pz = squark->momentum().z();
                  event_particles.push_back(Psquark);


                  // Fill event particle vector of chain final states ( including the quark)
                  MomentumVector Pquark;
                  Pquark.id = quark->id();
                  Pquark.m = quark->momentum().m();
                  Pquark.e = quark->momentum().e();
                  Pquark.px = quark->momentum().x();
                  Pquark.py = quark->momentum().y();
                  Pquark.pz = quark->momentum().z();
                  event_particles.push_back(Pquark);
                
                  MomentumVector Plepton1;
                  Plepton1.id = lepton1->id();
                  Plepton1.m  = lepton1->momentum().m();
                  Plepton1.e  = lepton1->momentum().e();
                  Plepton1.px = lepton1->momentum().x();
                  Plepton1.py = lepton1->momentum().y();
                  Plepton1.pz = lepton1->momentum().z();
                  event_particles.push_back(Plepton1);
                
                  MomentumVector Plepton2;
                  Plepton2.id = lepton2->id();
                  Plepton2.m  = lepton2->momentum().m();
                  Plepton2.e  = lepton2->momentum().e();
                  Plepton2.px = lepton2->momentum().x();
                  Plepton2.py = lepton2->momentum().y();
                  Plepton2.pz = lepton2->momentum().z();
                  event_particles.push_back(Plepton2);
                
                  MomentumVector Pneutralino1;
                  Pneutralino1.id = neutralino1->id();
                  Pneutralino1.m  = neutralino1->momentum().m();
                  Pneutralino1.e  = neutralino1->momentum().e();
                  Pneutralino1.px = neutralino1->momentum().x();
                  Pneutralino1.py = neutralino1->momentum().y();
                  Pneutralino1.pz = neutralino1->momentum().z();
                  event_particles.push_back(Pneutralino1);

                  // cout << "filled all momentum vectors" << endl;
                } // end if slepton has the right children
                else
                {
                  // cout << "slepton didn't have the right children" << endl;
                  chain = false;
                }
              } // end if getChildren2(slepton)
              else
              {
                // cout << "slepton didn't have children! (getChildren = false)" << endl;
                chain = false;
              }
            } // end if neutralino2 has the right children
            else
            {
              // cout << "neutralino 2 didn't have the right children" << endl;
              chain = false;
            }
          } // end if getChildren2(neutralino2)
          else
          {
            // cout << "neutralino2 didn't have children! (getChildren = false)" << endl;
            chain = false;
          }
        } // end if squark has the right children
        else
        {
          // cout << "squark didn't have the right children" << endl;
          chain = false;
        }
      } // end if getChildren2(squark)
      else
      {
        // cout << "squark didn't have children! (getChildren = false)" << endl;
        chain = false;
      }
    } // end if both outgoing primaries are squarks
  }
  else if ((abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) == 1000021) || (abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004 && abs(squarks[0]->id()) == 1000021))
  {
    // cout << "there weren't two squarks" << endl;
    gluino = true;
    chain = false;
  }
  else
  {
    // cout << "there weren't two squarks" << endl;
    chain = false;
  }



  if (chain)
  {
    numEvents++;
    if (momentumIsConserved)
    {
      numEventsMomentumConservation++;
    // Updated 20150116 to only write event to file if it passes momentum conservation tests

    // Write run number
    textOutput << "=== Chain event number " << numEvents << ", total event number " << numEventsTotal << " ===" << endl;
    // cout << "== Correct event number " << numEvents << ", total event number " << numEventsTotal << endl;
    //textOutput << "Number of final state particles: " << particles.size() << endl;


    // Write parseable 4-momenta to file
    textOutput << event_particles[0].id << "\t"
               << event_particles[0].px << "\t"
               << event_particles[0].py << "\t"
               << event_particles[0].pz << "\t"
               << event_particles[0].e  << "\t"
               << event_particles[0].m  << "\t"
               << endl;
    textOutput << event_particles[1].id << "\t"
               << event_particles[1].px << "\t"
               << event_particles[1].py << "\t"
               << event_particles[1].pz << "\t"
               << event_particles[1].e  << "\t"
               << event_particles[1].m  << "\t"
               << endl;
    textOutput << event_particles[2].id << "\t"
               << event_particles[2].px << "\t"
               << event_particles[2].py << "\t"
               << event_particles[2].pz << "\t"
               << event_particles[2].e  << "\t"
               << event_particles[2].m  << "\t"
               << endl;
    textOutput << event_particles[3].id << "\t"
               << event_particles[3].px << "\t"
               << event_particles[3].py << "\t"
               << event_particles[3].pz << "\t"
               << event_particles[3].e  << "\t"
               << event_particles[3].m  << "\t"
               << endl;
    textOutput << event_particles[4].id << "\t"
               << event_particles[4].px << "\t"
               << event_particles[4].py << "\t"
               << event_particles[4].pz << "\t"
               << event_particles[4].e  << "\t"
               << event_particles[4].m  << "\t"
               << endl;
    textOutput << event_particles[5].id << "\t"
               << event_particles[5].px << "\t"
               << event_particles[5].py << "\t"
               << event_particles[5].pz << "\t"
               << event_particles[5].e  << "\t"
               << event_particles[5].m  << "\t"
               << endl;
    textOutput << event_particles[6].id << "\t"
               << event_particles[6].px << "\t"
               << event_particles[6].py << "\t"
               << event_particles[6].pz << "\t"
               << event_particles[6].e  << "\t"
               << event_particles[6].m  << "\t"
               << endl;
    textOutput << event_particles[7].id << "\t"
               << event_particles[7].px << "\t"
               << event_particles[7].py << "\t"
               << event_particles[7].pz << "\t"
               << event_particles[7].e  << "\t"
               << event_particles[7].m  << "\t"
               << endl;
    // DEBUG: If we are writing the squark momenta to file, then the event_particles vector
    // contains two more entries which need to be written
    textOutput << event_particles[8].id << "\t"
               << event_particles[8].px << "\t"
               << event_particles[8].py << "\t"
               << event_particles[8].pz << "\t"
               << event_particles[8].e  << "\t"
               << event_particles[8].m  << "\t"
               << endl;
    textOutput << event_particles[9].id << "\t"
               << event_particles[9].px << "\t"
               << event_particles[9].py << "\t"
               << event_particles[9].pz << "\t"
               << event_particles[9].e  << "\t"
               << event_particles[9].m  << "\t"
               << endl;

    } // end if momentumIsConserved
  } // end if chain
  else if (gluino)
  {
   // textOutput << "== Total event number " << numEventsTotal << " (failed veto, gluino) ===" << endl;
  }
  else
  {
   // textOutput << "== Total event number " << numEventsTotal << " (failed veto) ===" << endl;
  }


  textOutput.close();  
}


LorentzRotation Foo::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void Foo::analyze(const tPVector & particles, double weight) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void Foo::analyze(tPPtr, double weight) {}

void Foo::dofinish() {
  AnalysisHandler::dofinish();
  // *** ATTENTION *** Normalize and post-process histograms here.
  ofstream textOutput;
  // ofstream adOut;
  textOutput.open("LHC-MSSM-analysis.log", ios::out | ios::app);
  //adOut.open("ad_data.dat", ios::out);
  textOutput << "Total number of events including wrong ones: " << numEventsTotal << endl;
  textOutput << "Number of events passing veto tests: " << numEvents << endl;
  textOutput << "Number of events passing veto and momentum conservation tests: " << numEventsMomentumConservation << endl;
  // textOutput << "Note that three-body decays of neutralino2 and slepton with gamma are accepted and gamma momentum discarded." << endl;
  textOutput << "Note that only events with full momentum conservation in the squark -> finals process is accepted" << endl;

  textOutput.close();  

  cout << "Total number of events including wrong ones: " << numEventsTotal << endl;
  cout << "Number of events passing veto tests: " << numEvents << endl;
  cout << "Number of events passing veto and momentum conservation tests: " << numEventsMomentumConservation << endl;
  // cout << "Note that three-body decays of neutralino2 and slepton with gamma are accepted and gamma momentum discarded." << endl;
  cout << "Note that only events with full momentum conservation in the squark -> finals process is accepted" << endl;
  // adOut.close();
}

void Foo::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
  numEvents = 0;
  numEventsTotal = 0;
  numEventsMomentumConservation = 0;
  ofstream textOutput;
  textOutput.open("LHC-MSSM-analysis.log", ios::trunc);
  textOutput.close();
}


IBPtr Foo::clone() const {
  return new_ptr(*this);
}

IBPtr Foo::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<Foo,AnalysisHandler>
  describeJorgenFoo("Jorgen::Foo", "Foo.so");

void Foo::Init() {

  static ClassDocumentation<Foo> documentation
    ("There is no documentation for the Foo class");

}

