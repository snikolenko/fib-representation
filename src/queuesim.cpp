//============================================================================
// Name        : queuesim.cpp
// Author      : Sergey Nikolenko
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/program_options.hpp>
#include "genpacket.hpp"
#include "queue.hpp"
#include "new_electric.hpp"

// #define DEBUG

using namespace std;
namespace po = boost::program_options;

double ELECTRIC_RADIUS = 0.5;


int num_runs = 5000000;
int max_threads = 1;

int min_usernum = 5;
int max_usernum = 25;
int default_usernum = 10;
int min_psnum = 5;
int max_psnum = 15;
int default_psnum = 10;
int min_maxlen = 3;
int max_maxlen = 100;
int default_maxlen = 10;
int max_capnum = 20;

int default_firstcap = 1;

string caida_infile = "";

// networking variables

int k_min = 10;
int k_max = 10;
int b_min = 10;
int b_max = 10;
int b_step = 1;
int c_min = 1;
int c_max = 1;
int l_min = 1;
int total_l = 1;


double default_beta = 1;
double beta = default_beta;


char mode = 'i';

bool only_two_inputs = false;
bool two_different_inputs = false;
bool no_smart = false;
bool no_depth_two = false;


void one_tick (int tick_num, PacketGenerator<int> * pg, vector<QueueContainer *> & v, uint k, bool add_packets = true) {
	D("\t\ttick");
	vector<IntPacket> p;
	if (add_packets) {
		p = pg->gen_packets();
		for (vector<IntPacket>::iterator it = p.begin(); it != p.end(); ++it) {
			it->setArrival(tick_num);
			D((*it) << " ");
		}
	}
	for (vector<QueueContainer *>::iterator it = v.begin(); it != v.end(); ++it) {
		D("\n\n\t### " << (*it)->type() << " ###");
		#ifdef DEBUG
			ostringstream s; for( vector<IntPacket>::const_iterator itt = p.begin(); itt != p.end(); ++itt ) s << *itt << " ";
			D("adding packets: " << s.str());
		#endif
		(*it)->add_packets( p );
		#ifdef DEBUG
			(*it)->print_queue("Added:");
		#endif
		(*it)->process(tick_num);
		#ifdef DEBUG
			(*it)->print_queue("Procd:");
		#endif
	}
}

void populate_manyqueue_vector ( vector<QueueContainer *> & v, uint k, uint min_B, uint max_B, uint B_step, uint min_C, uint max_C, double beta ) {
	v.clear();
	D("Populating vector, B=[" << min_B << ", " << max_B << "]");
	for (uint B = min_B; B <= max_B; B+=B_step ) {
		for (uint C = min_C; C <= max_C; ++C ) {
			D(k << "\t" << B << "\t" << C);
			v.push_back(new IntQueueContainer("PQ", k, B, C*k, false, 1));
			v.push_back(new SMQContainer(new NHSTSharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new NESTSharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new NHDTSharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new LQDSharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new BPDSharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new BPD2SharedMultiQueue(k, B, C)));
			v.push_back(new SMQContainer(new LWDSharedMultiQueue(k, B, C)));
		}
	}
}

void populate_vector ( vector<QueueContainer *> & v, uint k, uint min_B, uint max_B, uint B_step, uint min_C, uint max_C, double beta ) {
	v.clear();
	D("Populating vector, B=[" << min_B << ", " << max_B << "]");
	for (uint B = min_B; B <= max_B; B+=B_step ) {
		for (uint C = min_C; C <= max_C; ++C ) {
			D(k << "\t" << B << "\t" << C);
			if (beta > 3) {
				for (double bt = 1.0; bt < beta; bt += 0.1) {
					v.push_back(new IntQueueContainer("PQ", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("LPQ", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("FIFO", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("LFIFO", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("2LFIFO", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("RFIFO", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("LRFIFO", k, B, C, false, bt));
					v.push_back(new IntQueueContainer("RevPQ", k, B, C, false, bt));
				}
			} else {

			v.push_back(new IntQueueContainer("PQ", k, B, C, false));
			v.push_back(new IntQueueContainer("LPQ", k, B, C, false));
			v.push_back(new IntQueueContainer("NFIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("FIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("LFIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("2LFIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("RFIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("LRFIFO", k, B, C, false));
			v.push_back(new IntQueueContainer("RevPQ", k, B, C, false));

			if (beta != 1) {
				v.push_back(new IntQueueContainer("PQ", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("LPQ", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("FIFO", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("LFIFO", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("2LFIFO", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("RFIFO", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("LRFIFO", k, B, C, false, beta));
				v.push_back(new IntQueueContainer("RevPQ", k, B, C, false, beta));
			}
			}
		}
	}
}

void populate_lengthaware_vector ( vector<QueueContainer *> & v, uint k, uint min_B, uint max_B, uint B_step, uint min_C, uint max_C, uint min_L, uint max_L ) {
	v.clear();
	for (uint B = min_B; B <= max_B; B+=B_step ) {
		for (uint C = min_C; C <= max_C; ++C ) {
			for (uint L = min_L; L <= max_L; ++L ) {
				D("Generating: " << k << "\t" << B << "\t" << C << "\t" << L);
				v.push_back( new IntQueueContainer("FOPTUW", k, B, C, L, false) );
				v.push_back( new IntQueueContainer("POLen", k, B, C, L, false) );
				v.push_back( new IntQueueContainer("POWork", k, B, C, L, false) );
				v.push_back( new IntQueueContainer("POValue", k, B, C, L, false) );
			}
		}
	}
}

// packet generator parameters
double off_on_prob, on_off_prob, lmb, large_lmb;
int num_streams;
int max_length = 10;

void network_runsim(int k, int b_min, int b_max, int b_step, int c_min, int c_max, double beta, int num_runs) {
	vector<QueueContainer *> v;

	D("\toff_on_prob=" << off_on_prob << "\ton_off_prob=" << on_off_prob << "\tlmb=" << lmb << "\tlarge_lmb=" << large_lmb << "\tnum_streams=" << num_streams);
	// Globals::uniintK = boost::random::uniform_int_distribution<>(1, k);
	// Globals::uniintL = boost::random::uniform_int_distribution<>(1, 1);
	PacketGenerator<int> * pg = (PacketGenerator<int>*)(new MMPPVectorPoissonPacketGenerator<int>(k, off_on_prob, on_off_prob, lmb, large_lmb, num_streams));
	// PacketGenerator<int> * pg = (PacketGenerator<int>*)(new MMPPVectorPoissonPacketGenerator<int>(k, max_length, off_on_prob, on_off_prob, lmb, large_lmb, num_streams));
	D("\tgenerator created");
	populate_manyqueue_vector(v, k, b_min, b_max, b_step, c_min, c_max, beta);
	// populate_vector(v, k, b_min, b_max, b_step, c_min, c_max, beta);
	// populate_lengthaware_vector(v, k, b_min, b_max, b_step, c_min, c_max, max_length, max_length);
	for (int i=0; i < num_runs; ++i) {
		one_tick(i, pg, v, k);
		// emulate flush and restart
		if (i > 0 && i % 100000 == 0) {
			for (int j=0; j<k*b_max; ++j) {
				one_tick(i+j, pg, v, k, false);
			}
			i = i+k*b_max-1;
		}
	}
	#pragma omp critical
	{
		cout << "lambda=" << lmb << "\tlambda_on=" << large_lmb << "\tk=" << k <<
				"\ttotal_packets=" << pg->total_packets <<
				"\ttotal_packets_len=" << pg->total_packets_len <<
				"\tnum_streams=" << num_streams << endl;
		double max_total_length = v[0]->final_total_processed_length();
		for (vector<QueueContainer *>::iterator it = v.begin(); it != v.end(); ++it) {
			cout << (*it)->type().c_str() << "\tlambda=" << large_lmb << "\tk=" << (*it)->k() << "\tB=" << (*it)->B()
					<< "\tC=" << (*it)->C() << "\tbeta=" << std::setprecision(2) << (*it)->beta() << "\t"
					<< (*it)->final_total_processed() << "\t"
					<< int((*it)->final_total_processed_length()) << "\t"
					// << std::setprecision(5) << double((*it)->final_total_processed_length() / max_total_length) << "\t"
					<< (*it)->total_admitted() << "\t"
					<< std::setprecision(5) << ((*it)->total_delay() / (double)(*it)->total_processed()) << "\t"
					<< std::setprecision(5) << sqrt((*it)->total_squared_delay()) / (double)(*it)->total_processed()
					<< endl;
		}
	}
}


using namespace boost;

int network_main() {
	srand((unsigned)time(0));

	double start_l = 0;
	double step_l = 0;
	// double large_l_start = 5;
	// double large_l_step = 10;
	double large_l_start = 0.025;
	double large_l_step = 0.01;
	double default_lmb = 0.25;

	num_streams = 500;
	off_on_prob = 0.01;
	on_off_prob = 0.2;
	max_length = default_maxlen;

	size_t SIM_SIZE = num_runs;

	vector< vector< vector<QueueContainer *> > > v(total_l);
	vector< vector<int> > total_packets(total_l);
	vector< vector<int> > total_packets_len(total_l);

	omp_set_num_threads(2);

	for ( int l=0; l < total_l; ++l ) {
		lmb = start_l + l * step_l;
		large_lmb = large_l_start + l * large_l_step;
		if (total_l == 1) large_lmb = default_lmb;
		for ( int k=k_min; k <= k_max; k += 2 ) {
			// num_streams = k * max_length;
			// large_lmb = (large_l_start + l * large_l_step) / (double)num_streams;
			D("Experiment: k=" << k <<"\tb_min=" << b_min << "\tb_max=" << b_max << "\tb_step=" << b_step << "\tbeta=" << beta << "\truns=" << SIM_SIZE);
			network_runsim(k, b_min, b_max, b_step, c_min, c_max, beta, SIM_SIZE);
		}
	}
	return 0;
}

void electric_populate( vector<Algorithm *> & v ) {
	v.clear();
	if (only_two_inputs) {
		v.push_back(new SLDAlgorithm("SLD"));
		v.push_back(new SLPAlgorithm("SLP"));
		v.push_back(new SBRAlgorithm("SBP"));
		for (size_t i=0; i<v.size(); ++i) v[i]->shared_port = true;
	} else {
		v.push_back(new LDAlgorithm("LD"));
		v.push_back(new RDAlgorithm("RD"));
		v.push_back(new SOPAlgorithm("SOP"));
		v.push_back(new LSOPAlgorithm("LSOP"));
	}
}

void electric_one_try(size_t userNum, size_t psNum, int maxlen, vector<Algorithm *> & v, bool limit_two, int first_port_cap = 1, double radius = ELECTRIC_RADIUS) {
	vector<User> users;
	vector<PowerStation> pstations;
	for (size_t i=0; i<psNum; ++i) {
		pstations.push_back(PowerStation());
	}
	for (size_t i=0; i<userNum; ++i) {
		users.push_back(User(maxlen, pstations, radius, limit_two, two_different_inputs));
	}

	for (size_t i = 0; i < v.size(); ++i) {
		v[i]->process(users, first_port_cap);
	}
}

void electric_runsim(size_t userNum, size_t psNum, int maxlen, int num_runs, bool limit_two, int first_port_cap = 1, double radius = ELECTRIC_RADIUS) {
	vector<Algorithm *> v;
	electric_populate(v);

	for (int i=0; i < num_runs; ++i) {
		electric_one_try(userNum, psNum, maxlen, v, limit_two, first_port_cap, radius);
	}

	#pragma omp critical
	{
		cout << "usernum=" << userNum << "\tpsnum=" << psNum << "\tmaxlen=" << maxlen << "\tfirstcap=" << first_port_cap << "\tradius=" << radius << endl;
		for (size_t i = 0; i < v.size(); ++i) {
			v[i]->print_stats();
		}
	}
}

int electric_main() {
	srand((unsigned)time(NULL));
	omp_set_num_threads(max_threads);

	if (mode == 'i') {
		#pragma omp parallel for
		for (int psnum = min_psnum; psnum <= max_psnum; psnum += 5) {
			electric_runsim(default_usernum, psnum, default_maxlen, num_runs, only_two_inputs, default_firstcap);
		}
	}

	if (mode == 'd') {
		#pragma omp parallel for
		for (int usnum = min_usernum; usnum <= max_usernum; usnum += 25) {
			electric_runsim(usnum, default_psnum, default_maxlen, num_runs, only_two_inputs, default_firstcap);
		}
	}

	if (mode == 'c') {
		#pragma omp parallel for
		for (int capnum = 1; capnum <= max_capnum; capnum += 1) {
			electric_runsim(default_usernum, default_psnum, default_maxlen, num_runs, only_two_inputs, capnum);
		}
	}

	if (mode == 'l') {
		#pragma omp parallel for
		for (int mlen = min_maxlen; mlen <= max_maxlen; mlen += 50) {
			electric_runsim(default_usernum, default_psnum, mlen, num_runs, only_two_inputs, default_firstcap);
		}
	}

	if (mode == 'r') {
		#pragma omp parallel for
		for (int radius = 1; radius < 20; radius++) {
			electric_runsim(default_usernum, default_psnum, default_maxlen, num_runs, only_two_inputs, default_firstcap, radius * 0.05);
		}
	}

	return 0;
}

int main(int argc, char* argv[]) {
	//return network_main();
	// Declare the supported options.

	srand(time(NULL));

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,?", "produce help message")
	    ("max_threads,t", po::value<int>(&max_threads)->default_value(1), "number of threads to run")
	    ("mode,m", po::value<char>(&mode)->default_value('i'), "type of experiment ('i'nput, 'd'emand, 'l'ength for electric; k, b, c, be[t]a for network)")
	    ("num_runs", po::value<int>(&num_runs)->default_value(10000), "number of simulation runs")
	    ("only_two,2", "two inputs per user max")
	    ("two_different", "two per user with different inputs")
	    ("no_smart,n", "no smart algorithms (much faster)")
	    ("no_dtwo", "no depth two algorithms (even faster)")
	    ("min_usernum", po::value<int>(&min_usernum)->default_value(5), "min number of users")
	    ("max_usernum", po::value<int>(&max_usernum)->default_value(25), "max number of users")
	    ("usernum", po::value<int>(&default_usernum)->default_value(10), "default number of users (for experiments w.r.t. other vars)")
	    ("min_psnum", po::value<int>(&min_psnum)->default_value(5), "min number of inputs")
	    ("max_psnum", po::value<int>(&max_psnum)->default_value(15), "max number of inputs")
	    ("psnum", po::value<int>(&default_psnum)->default_value(10), "default number of ports (for experiments w.r.t. other vars)")
	    ("min_maxlen", po::value<int>(&min_maxlen)->default_value(3), "min maximal demand length")
	    ("max_maxlen", po::value<int>(&max_maxlen)->default_value(100), "max maximal demand length")
	    ("maxlen", po::value<int>(&default_maxlen)->default_value(10), "default maximal demand length (for experiments w.r.t. other vars)")
	    ("firstcap", po::value<int>(&default_firstcap)->default_value(1), "default first port capacity")
	    ("k_min", po::value<int>(&k_min)->default_value(10), "min k")
	    ("k_max", po::value<int>(&k_max)->default_value(10), "max k")
	    ("b_min", po::value<int>(&b_min)->default_value(50), "min b")
	    ("b_max", po::value<int>(&b_max)->default_value(50), "max b")
	    ("b_step", po::value<int>(&b_step)->default_value(1), "step b")
	    ("c_min", po::value<int>(&c_min)->default_value(1), "min c")
	    ("c_max", po::value<int>(&c_max)->default_value(1), "max c")
	    ("l_min", po::value<int>(&l_min)->default_value(1), "min l")
	    ("total_l", po::value<int>(&total_l)->default_value(1), "min l")
	    ("beta", po::value<double>(&beta)->default_value(1.0), "beta")
	    ("caida", po::value<string>(&caida_infile)->default_value(""), "caida traces input")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	if (vm.count("only_two")) only_two_inputs = true;
	if (vm.count("no_smart")) no_smart = true;
	if (vm.count("no_depth_two")) no_depth_two = true;

	if (caida_infile != "") {
		cout << "Will use CAIDA traces from " << caida_infile << endl;
	}

	return network_main();
	// return electric_main();
}

