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
using namespace arma;

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

PPtr getFinalInstance(PPtr p)
{
  //int id = p->id();
  ParticleVector children = p->children();
  if(children.size()==0) 
  {
    return p;
  }
  else if (children[0]->id() != p->id())
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
  if (children.size()>1)
  {
    return p;
  }
  else if (children.size()==0)
  {
    cout << "final state quark! id = " << p->id() << endl;
    return p;
  }
  else
  {
    return getOffshellQuark(children[0]);
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

  /** Get the final-state particles */
  tPVector particles=event->getFinalState();

  /** Make and open text output file */
  ofstream textOutput;
  // ofstream adOut;
  textOutput.open("LHC-MSSM-analysis.log", ios::out | ios::app);

  // Declare variables 
  vector<MomentumVector> event_particles; // Save all chain particles
  bool chain = true; // Assume that both chains are present until proven wrong by if-tests
  // Get primary subprocess (squarks):
  tSubProPtr primarysubprocess = event->primarySubProcess();
  ParticleVector squarks = primarysubprocess->outgoing();

  // Check that we produced two left-handed squarks of first
  // or second generation
  if (abs(squarks[0]->id()) >= 1000001 && abs(squarks[0]->id()) <= 1000004 && abs(squarks[1]->id()) >= 1000001 && abs(squarks[1]->id()) <= 1000004)
  {
    cout << "loop over squarks" << endl;
    // Loop over squarks
    for (ParticleVector::const_iterator squarkiterator = squarks.begin(); squarkiterator != squarks.end(); ++squarkiterator)
    {
      cout << "entered squark loop" << endl;
      // Get squark children -- we know there have to be some, since it's a squark
      ParticleVector squarkchildren;
      if (getChildren2(*squarkiterator, squarkchildren))
      {
        cout << "getchildren2(squarks)" << endl;
        if ( squarkchildren.size()==2 && ( (abs(squarkchildren[0]->id()) >=1 && abs(squarkchildren[0]->id()) <=4 && squarkchildren[1]->id() == 1000023) || (abs(squarkchildren[1]->id()) >=1 && abs(squarkchildren[1]->id()) <=4 && squarkchildren[0]->id() == 1000023) ) )
        {
          cout << "squark has the right children" << endl;
          // Check which child is which, make sure the quark is showered
          PPtr neutralino2, quark;
          if (abs(squarkchildren[0]->id()) > abs(squarkchildren[1]->id()))
          {
          //  cout << "trying to get offshell quark of id " << squarkchildren[1]->id() << endl;
            neutralino2 = squarkchildren[0];
            quark = getFinalInstance(squarkchildren[1]);
          }
          else
          {
           // cout << "trying to get offshell quark of id " << squarkchildren[0]->id() << endl;
            neutralino2 = squarkchildren[1];
            quark = getFinalInstance(squarkchildren[0]);
          }

          cout << "got both squark children" << endl;
          // Get neutralino2 children -- has to have children since we know it's a neutralino2
          ParticleVector neutralino2children;
          if (getChildren2(neutralino2,neutralino2children))
          {
            cout << "getchildren2(neutralino)" << endl;
            if ( 
              neutralino2children.size()==2 && ( ( (abs(neutralino2children[0]->id()) ==11 || abs(neutralino2children[0]->id()) ==13) && ( abs(neutralino2children[1]->id()) == 1000011 || abs(neutralino2children[1]->id()) == 1000013 || abs(neutralino2children[1]->id()) == 2000011 || abs(neutralino2children[1]->id()) == 2000013 ) ) || ( (abs(neutralino2children[1]->id()) ==11 || abs(neutralino2children[1]->id()) ==13) && ( abs(neutralino2children[0]->id()) == 1000011 || abs(neutralino2children[0]->id()) == 1000013 || abs(neutralino2children[0]->id()) == 2000011 || abs(neutralino2children[0]->id()) == 2000013 ) ) ) 
              )
            {
              cout << "neutralino2 has the right children" << endl;
              // Check which child is which
              PPtr slepton, lepton1;
              if (abs(neutralino2children[0]->id()) > abs(neutralino2children[1]->id()))
              {
                slepton = neutralino2children[0];
                lepton1 = neutralino2children[1];
              }
              else
              {
                slepton = neutralino2children[1];
                lepton1 = neutralino2children[0];      
              }

              cout << "got both neutralino2 children" << endl;
              // Get slepton children -- has to have children since we know it's a slepton
              ParticleVector sleptonchildren;
              if (getChildren2(slepton, sleptonchildren))
              {
                cout << "getChildren2(slepton)" << endl;
                if ( 
                  sleptonchildren.size()==2 && ( ( (abs(sleptonchildren[0]->id()) ==11 || abs(sleptonchildren[0]->id()) ==13) && abs(sleptonchildren[1]->id()) == 1000022 ) || ( (abs(sleptonchildren[1]->id()) ==11 || abs(sleptonchildren[1]->id()) ==13) && abs(sleptonchildren[0]->id()) == 1000022 ) ) 
                  )
                {
                  cout << "slepton has the right children" << endl;
                  PPtr neutralino1, lepton2;
                  if (abs(sleptonchildren[0]->id()) > abs(sleptonchildren[1]->id()))
                  {
                    neutralino1 = sleptonchildren[0];
                    lepton2 = sleptonchildren[1];
                  }
                  else
                  {
                    neutralino1 = sleptonchildren[1];
                    lepton2 = sleptonchildren[0];
                  }

                  cout << "got both slepton children" << endl;



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

                  cout << "filled all momentum vectors" << endl;
                } // end if slepton has the right children
                else
                {
                  cout << "slepton didn't have the right children" << endl;
                  chain = false;
                }
              } // end if getChildren2(slepton)
              else
              {
                cout << "slepton didn't have children! (getChildren = false)" << endl;
                chain = false;
              }
            } // end if neutralino2 has the right children
            else
            {
              cout << "neutralino 2 didn't have the right children" << endl;
              chain = false;
            }
          } // end if getChildren2(neutralino2)
          else
          {
            cout << "neutralino2 didn't have children! (getChildren = false)" << endl;
            chain = false;
          }
        } // end if squark has the right children
        else
        {
          cout << "squark didn't have the right children" << endl;
          chain = false;
        }
      } // end if getChildren2(squark)
      else
      {
        cout << "squark didn't have children! (getChildren = false)" << endl;
        chain = false;
      }
    } // end if both outgoing primaries are squarks
  }
  else
  {
    cout << "there weren't two squarks" << endl;
    chain = false;
  }

  if (chain == true)
  {
    numEvents++;
    // Write run number
    textOutput << "=== Chain event number " << numEvents << ", total event number " << numEventsTotal << " ===" << endl;
    cout << "== Correct event number " << numEvents << ", total event number " << numEventsTotal << endl;
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
  }
  else
  {
    cout << "== Total event number " << numEventsTotal << endl;
  }


  textOutput.close();  
  // adOut.close();
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

  textOutput.close();  
  // adOut.close();
}

void Foo::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
  numEvents = 0;
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

