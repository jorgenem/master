// PART OF PROGRAM: Method to read events from file


struct MomentumVector
{
  int id;
  vector<double> p;
};


void best_fit( ... )
{


	int N = Nbins*Nevents;
	vector<bool> all_leptons_equal_list;
	
	string line;
	// ifstream events ("../events/on-shell_decay_squarks_at_rest_10000_events.dat");
	ifstream events ("../events/Pythia_cascade_events_no_ISR_or_FSR_20150120_only_opposite_flavour_leptons.dat");
	// ifstream events ("../events/Pythia_cascade_10000_events_everything_turned_on_20150210_only_opposite_flavour_leptons.dat");
	// ifstream events ("../events/Herwig_chain_20150116_with_gluinos_and_no_threebody_decay_and_discarded_momentum-nonconservation_GeV-corrected_only_opposite_flavour_leptons.dat");
	
	if (events.is_open())
	{
	
		for (int iEvent = 0; iEvent < N; iEvent++) 
		{
			vector<string> particle;
			MomentumVector p1, p2, p3, p4, p5, p6, p7, p8;
	
	
			for (int iParticle = 0; iParticle < 9; iParticle++)
			{
				getline(events,line);
				istringstream iss(line);
				particle = {istream_iterator<string>{iss}, istream_iterator<string>{}};
	
				if (iParticle == 1)
				{	
					p1.id = stoi(particle[0]);
					p1.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 2)
				{
					p2.id = stoi(particle[0]);
					p2.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 3)
				{
					p3.id = stoi(particle[0]);
					p3.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 4)
				{
					p4.id = stoi(particle[0]);
					p4.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 5)
				{
					p5.id = stoi(particle[0]);
					p5.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 6)
				{
					p6.id = stoi(particle[0]);
					p6.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 7)
				{
					p7.id = stoi(particle[0]);
					p7.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
				if (iParticle == 8)
				{
					p8.id = stoi(particle[0]);
					p8.p = {stod(particle[1]), stod(particle[2]), stod(particle[3]), stod(particle[4])};
				}
			}
	
			all_leptons_equal_list.push_back((abs(p2.id)==abs(p6.id)));
	
	
	
			// END FOR loop over events
		}
	    events.close();
	}
	else cout << "Unable to open file"; 	

}