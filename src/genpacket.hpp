/*
 * genpacket.hpp
 *
 *  Created on: Jun 8, 2012
 *      Author: snikolenko
 */

#ifndef GENPACKET_HPP_
#define GENPACKET_HPP_

#include <vector>
#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "packet.hpp"

using namespace std;

bool bernoulli( double prob );

template <typename T> class PacketGenerator {
public:
	PacketGenerator(int ik) : k(ik), L(1), total_packets(0), total_packets_len(0) {};
	PacketGenerator(int ik, int iL) : k(ik), L(iL), total_packets(0), total_packets_len(0) {};

	vector<Packet<T> > gen_packets() {
		vector<Packet<T> > res = internal_gen_packets();
		total_packets += res.size();
		for (typename vector<Packet<T> >::const_iterator it = res.begin(); it != res.end(); ++it) {
			total_packets_len += it->l;
		}
		return res;
	}

	virtual vector<Packet<T> > internal_gen_packets() {
		int nump = gen_n_packets();
		return gen_several_packets(nump);
	}

	virtual Packet<T> gen_packet() {
		return Packet<T>( rand() % k + 1, rand() % L + 1 );
	}

	virtual vector<Packet<T> > gen_several_packets(int nump) {
		vector<Packet<T> > res;
		for (int i=0; i < nump; ++i) {
			res.push_back( gen_packet() );
		}
		return res;
	}

	virtual int gen_n_packets() {
		return 0;
	}

	int k;
	int L;
	int total_packets;
	int total_packets_len;
};

template <typename T> class MMPPUniformPacketGenerator : public PacketGenerator<T> {
public:
	MMPPUniformPacketGenerator<T>( int ik, double swprob, double loff, int minon, int maxon ) : PacketGenerator<T>(ik),
			on(false), switch_prob(swprob), lambda_off(loff), min_on(minon), max_on(maxon) {};
	MMPPUniformPacketGenerator<T>( int ik, int iL, double swprob, double loff, int minon, int maxon ) : PacketGenerator<T>(ik, iL),
			on(false), switch_prob(swprob), lambda_off(loff), min_on(minon), max_on(maxon) {};

	vector<Packet<T> > internal_gen_packets();

	bool on;
	double switch_prob;
	double lambda_off;
	int min_on;
	int max_on;
};

template <typename T> class UniformPacketGenerator : public PacketGenerator<T> {
public:
	UniformPacketGenerator<T>( int ik, int minnum, int maxnum ) : PacketGenerator<T>(ik),
			min_num(minnum), max_num(maxnum) {};
	UniformPacketGenerator<T>( int ik, int iL, int minnum, int maxnum ) : PacketGenerator<T>(ik, iL),
			min_num(minnum), max_num(maxnum) {};

	int gen_n_packets();

	int min_num;
	int max_num;
};

template <typename T> class PoissonPacketGenerator : public PacketGenerator<T> {
public:
	PoissonPacketGenerator<T>( int ik, double ilambda) : PacketGenerator<T>(ik),
			lambda(ilambda), gen((unsigned)time(0)), pois(ilambda > 0 ? ilambda : 1), rpois(gen, pois) { };
	PoissonPacketGenerator<T>( int ik, int iL, double ilambda) : PacketGenerator<T>(ik, iL),
			lambda(ilambda), gen((unsigned)time(0)), pois(ilambda > 0 ? ilambda : 1), rpois(gen, pois) { };

	int gen_n_packets();

	double lambda;
	boost::mt19937 gen;
	boost::poisson_distribution<int, double> pois;
	boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int, double> > rpois;
};

int get_random_int(int min, int max);

template <typename T> class MMPPoissonPacketGenerator : public PacketGenerator<T> {
public:
	MMPPoissonPacketGenerator<T>( int ik, double swprob, double swbprob, double loff, double lon) : PacketGenerator<T>(ik),
			k(ik), L(1), on(false), switch_prob(swprob), switch_back_prob(swbprob), pg_on(ik, lon), pg_off(ik, loff), gen(42u) {};
	MMPPoissonPacketGenerator<T>( int ik, int iL, double swprob, double swbprob, double loff, double lon) : PacketGenerator<T>(ik, iL),
			k(ik), L(1), on(false), switch_prob(swprob), switch_back_prob(swbprob), pg_on(ik, lon), pg_off(ik, loff), gen(42u) {
				// cout << "\tgenerated MMPP\n";
			};

	int gen_n_packets();

	virtual Packet<T> gen_packet() {
		int kk = (k == 1) ? 1 : get_random_int(1, k);
		// cout << "\n\tgenerated " << kk << "\n";
		int ll = (L == 1) ? 1 : get_random_int(1, L);
		Packet<T> p = Packet<T>( kk, ll );
		return p;
	}

	int k, L;
	bool on;
	double switch_prob;
	double switch_back_prob;
	PoissonPacketGenerator<T> pg_on;
	PoissonPacketGenerator<T> pg_off;
	boost::mt19937 gen;
};


template <typename T> class MMPPVectorPoissonPacketGenerator : public PacketGenerator<T> {
public:
	MMPPVectorPoissonPacketGenerator<T>( int ik, double swprob, double swbprob, double loff, double lon, int n ) : PacketGenerator<T>(ik),
			on(false), switch_prob(swprob), switch_back_prob(swbprob), lambda_off(loff), lambda_on(lon), nstreams(n) {
		initVector();
	};
	MMPPVectorPoissonPacketGenerator <T>( int ik, int iL, double swprob, double swbprob, double loff, double lon, int n ) : PacketGenerator<T>(ik, iL),
			on(false), switch_prob(swprob), switch_back_prob(swbprob), lambda_off(loff), lambda_on(lon), nstreams(n) {
		initVector();
	};

	void initVector() {
		for (int i=0; i < nstreams; ++i) {
			// int cur_k = i % this->k + 1;
			// int cur_L = (i / this->k) % this->L + 1;
			// cout << "Stream " << i << ":  k=" << cur_k << "\tL=" << cur_L << endl;
			// v.push_back(new MMPPoissonPacketGenerator<T>(cur_k, cur_L, switch_prob, switch_back_prob, lambda_off, lambda_on));
			v.push_back(new MMPPoissonPacketGenerator<T>(this->k, this->L, switch_prob, switch_back_prob, lambda_off, lambda_on));
		}
	}

	int gen_n_packets() { return 0; }

	vector<Packet<T> > internal_gen_packets();

	vector<MMPPoissonPacketGenerator<T> *> v;
	bool on;
	double switch_prob;
	double switch_back_prob;
	double lambda_off;
	double lambda_on;
	int nstreams;
};


#endif /* GENPACKET_HPP_ */
