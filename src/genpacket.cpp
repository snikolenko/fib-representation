/*
 * genpacket.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: snikolenko
 */
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "genpacket.hpp"
#include "queue.hpp"
#include <stdlib.h>
#include <time.h>

bool bernoulli( double prob ) {
	return ( ( rand() / (double)RAND_MAX ) < prob );
}

int get_random_int(int min, int max) {
	return min + (rand() % (max-min+1));
}

template <typename T> vector<Packet<T> > MMPPUniformPacketGenerator<T>::internal_gen_packets() {
	if ( bernoulli(switch_prob) ) on = !on;
	D("gen_packets " << on);
	vector<Packet<T> > res;
	if (on) {
		int nump = rand() % (max_on - min_on) + min_on;
		D("\t" << nump << "packets");
		for (int i=0; i < nump; ++i) {
			res.push_back( this->gen_packet() );
		}
	} else {
		if ( lambda_off > 0 && bernoulli(lambda_off) ) {
			D("\t" << "1 packet in off state");
			res.push_back( this->gen_packet() );
		}
	}
	return res;
}

template <typename T> int UniformPacketGenerator<T>::gen_n_packets() {
	return ( rand() % (max_num - min_num) + min_num );
}

template <typename T> int PoissonPacketGenerator<T>::gen_n_packets() {
	return (lambda > 0 ? rpois() : 0);
}

template <typename T> int MMPPoissonPacketGenerator<T>::gen_n_packets() {
	if ( on && bernoulli(switch_back_prob) ) on = false;
	if ( !on && bernoulli(switch_prob) ) on = true;
	if (on) return pg_on.gen_n_packets();
	else return pg_off.gen_n_packets();
}

template <typename T, typename L> vector<Packet<T> > MMPPVectorPoissonPacketGenerator<T, L>::internal_gen_packets() {
	vector<Packet<T> > res;
	for (typename vector<L *>::iterator it = v.begin(); it != v.end(); ++it) {
		vector<Packet<T> > current_res = (*it)->internal_gen_packets();
		res.insert(res.end(), current_res.begin(), current_res.end());
	}
	//cout << " " << res.size();
	return res;
}

template class MMPPoissonPacketGenerator<int>;
template class MMPPVectorPoissonPacketGenerator<int, MMPPoissonPacketGenerator<int> >;
template class MMPPVectorPoissonPacketGenerator<int, MMPPoissonBiasedPacketGenerator<int> >;
template class MMPPVectorPoissonPacketGenerator<int, MMPPoissonTwoValuedBiasedPacketGenerator<int> >;
template class MMPPVectorPoissonPacketGenerator<int, MMPPoissonTwoValuedUniformPacketGenerator<int> >;
template class PoissonPacketGenerator<int>;
template class PacketGenerator<int>;
