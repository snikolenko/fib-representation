/*
 * queue.cpp
 *
 *  Created on: Jun 8, 2012
 *      Author: snikolenko
 */
#include "queue.hpp"

#include <iostream>
#include <sstream>
#include <algorithm>

template <typename T> void Queue<T>::init() {
	if (multiqueue) {
		for (size_t i=0; i < k; ++i ) {
			mq.push_back( Queue<T>("FIFO", k, B, 1, false) );
		}
		curq = 0;
	}
	preempting = true;
	if (type == "NFIFO") {
		preempting = false;
	}
	if (type == "2LFIFO") {
		leavepartiallyprocessed = true;
	}
	if (type == "PQ" || type == "LPQ" || type == "RevPQ") {
		sorting = true;
	}
	if (type == "LPO" || type == "LPQ" || type == "LFIFO" || type == "2LFIFO" || type == "LRFIFO") {
		lazy = true;
	}
	if (type == "RFIFO" || type == "LRFIFO") {
		recycling = true;
	}
	if (type == "RevPQ" || type == "RevFIFO") {
		reversing = true;
	}

	// length-aware policies
	if ( (type == "POLen") || (type == "POWork") || (type == "POValue") ) {
		preempting = true;
		sorting = true;
		lengthaware = true;
	}
	if ( (type == "FOPTUL") || (type == "FOPTUW") ) {
		preempting = true;
		sorting = true;
		lengthaware = true;
		fractional = true;
	}
}

template <typename T> void Queue<T>::increment_counters(int length, int delay, int squared_delay) {
	total_processed++;
	total_processed_length += length;
	total_delay += delay;
	total_squared_delay += squared_delay;
}

template <typename T> void Queue<T>::remove_head_packet(int tick_num) {
	typename vector< Packet<T> >::iterator it = this->reversing ? (q.end()-1) : q.begin();
	D("[" << tick_num << "]: processed " << (*it) << " adding " << (tick_num - it->arrival - (int)it->r + 1));
	increment_counters(it->l, tick_num - it->arrival - (int)it->r + 1, (tick_num - it->arrival) * (tick_num - it->arrival));
	q.erase(it);
	if ( type != "FOPTUW" ) {
		D(type << " processed a packet!\ntotal_proc=" << total_processed << "\ttotal_len=" << total_processed_length << "\ttotal_delay=" << total_delay << "\tlen_left=" << totallength());
	}
}

template <typename T> void Queue<T>::process_head_packet(int tick_num) {
	if ( q.empty() ) return;
	else {
		if (!fractional) {
			if (lazy) {
				// lazy packet processing
				if ( q[0].r != 1 && sendingout ) sendingout = false;
				if (sendingout) this->remove_head_packet(tick_num);
				else {
					size_t non_one = 0;
					for ( ; (non_one < q.size()) && (q[non_one] == 1); ++non_one );
					if (non_one == q.size()) {
						sendingout = true;
						this->remove_head_packet(tick_num);
					} else {
						q[non_one].dec();
					}
				}
			} else {
				size_t holindex = this->reversing ? (q.size()-1) : 0;
				// common integer packet processing
				if ( q[holindex] == 1 ) this->remove_head_packet(tick_num);
				else {
					q[holindex].dec();
					if ( recycling ) {
						this->recycle();
					}
				}
			}
		} else {
			// fractional packet processing
			float work_left = 1;
			while ( !q.empty() && work_left >= q[0].r ) {
				work_left -= q[0].r;
				this->remove_head_packet(tick_num);
			}
			if ( !q.empty() ) {
				q[0].r -= work_left;
			}
		}
	}
}

template <typename T> void Queue<T>::print_queue(string pref = "") {
	if (multiqueue) {
		D("Multiqueue state for " << type);
		for (size_t i=0; i<k; ++i) {
			mq[i].print_queue("   ");
		}
	} else {
		ostringstream s; s << pref << "\t";
		for (typename vector<Packet<T> >::const_iterator it = q.begin(); it != q.end(); ++it) {
			s << *it << " ";
		}
		D(s.str());
	}
}

template <typename T> void Queue<T>::process_subqueue(int tick_num, int i) {
	D("  " << type << ": processing subqueue " << (i+1));
	curq = i;
	T prevtotal = mq[i].total_processed;
	uint prevtotallen = mq[i].total_processed_length;
	uint prevtotaldelay = mq[i].total_delay;
	uint prevtotalsqdelay = mq[i].total_squared_delay;
	mq[i].process(tick_num);
	if ( mq[i].total_processed > prevtotal) {
		total_processed++;
		total_processed_length += mq[i].total_processed_length - prevtotallen;
		total_delay += mq[i].total_delay - prevtotaldelay;
		total_squared_delay += mq[i].total_squared_delay - prevtotalsqdelay;
		D(type << " processed a packet!\ntotal_proc=" << total_processed << "total_len=" << total_processed_length);
	}
}

template <typename T> void Queue<T>::process(int tick_num) {
	for (size_t j=0; j<C; ++j) {
		if (multiqueue) {
			process_multiqueue(tick_num);
			continue;
		}
		if ( q.empty() ) return;
		process_head_packet(tick_num);
	}
}

template <typename T> void Queue<T>::recycle() {
	if (q.size() > 1) {
		Packet<T> hol = q[0];
		q[0] = q[q.size()-1];
		q[q.size()-1] = hol;
	}
}

template <typename T> int Queue<T>::get_next_queue(size_t cur) {
	D("get_next_queue cur=" << cur);
	size_t first_i = (cur == k-1) ? 0 : cur+1;
	size_t i = first_i;
	while (true) {
		if (i == k) i=0;
		if ( mq[i].has_packets() ) return i;
		if (i == cur) break;
		++i;
	}
	D("		return " << 0);
	return 0;
}

template <typename T> void Queue<T>::process_multiqueue(int tick_num) {
	//TODO: pm
	int res = -1;
	if (type == "MQF") {
		res = 0;
		for ( size_t i=0; i < k; ++i ) {
			if (mq[i].has_packets() ) {
				res = i; break;
			}
		}
	} else if (type == "MaxQF") {
		res = k-1;
		for ( int i=(int)k-1; i >= 0; --i ) {
			if (mq[i].has_packets() ) {
				res = i; break;
			}
		}
	} else if (type == "CRR") {
		res = get_next_queue(curq);
	} else if (type == "PRR") {
		if (mq[curq].has_packets() && mq[curq].hol() < (curq+1) ) res = curq;
		else res = get_next_queue(curq);
		D("  PRR curq=" << curq+1 << "  res=" << res+1);
	} else if (type == "LQF") {
		res = max_element( mq.begin(), mq.end(), compare_queues_by_length ) - mq.begin();
		D("  LQF longest=" << res+1);
	} else if (type == "SQF") {
		size_t minind = 0;
		size_t minlen = (mq[0].has_packets() ? mq[0].size() : (B+1));
		for ( size_t i=1; i < k; ++i ) {
			if (mq[i].has_packets() && mq[i].size() < minlen ) {
				minlen = mq[i].size(); minind = i;
			}
		}
		res = minlen > 0 ? minind : 0;
		D("  SQF shortest=" << res+1);
	} else {
		res = 0;
	}
	D("Processing queue " << (res+1));
	process_subqueue(tick_num, res);
}

template <typename T> void Queue<T>::add_packets( const vector<IntPacket> & v ) {
	D("Adding " << v.size() << " packets");
	for ( vector<IntPacket>::const_iterator it = v.begin(); it != v.end(); ++it ) {
		add_packet( *it, false );
	}
	if (sorting) {
		doSort();
	}
	if (lengthaware && preempting) {
		#ifdef DEBUG
			print_queue("Before preempting: ");
		#endif
		int preempt_to = B-2*L+1;
		if ( (type == "FOPTUL") || (type == "FOPTUW") ) preempt_to = B;
		D("Preempting to " << preempt_to << "\tB=" << B << " L=" << L);
		// preemption for length-aware policies
		double len = totallength();
		while ( len > preempt_to + EPSILON ) {
			len -= q[q.size()-1].l;
			q.erase(q.begin()+(q.size()-1));
		}
	}
}

template <typename T> void Queue<T>::doSort() {
	if (!lengthaware) {
		sort(q.begin(), q.end());
	} else {
		if ( (type == "POLen") || (type == "FOPTUW") ) {
			sort(q.begin(), q.end(), sortLength<T>);
		}
		else if ((type == "POWork") || (type == "FOPTUL")) {
			sort(q.begin(), q.end(), sortWork<T>);
		}
		else if (type == "POValue") {
			sort(q.begin(), q.end(), sortValue<T>);
		}
	}
}

template <typename T> typename vector<Packet<T> >::iterator Queue<T>::get_max_to_preempt() {
	if (this->leavepartiallyprocessed) {
		typename vector<Packet<T> >::iterator res = q.end();
		for (typename vector<Packet<T> >::iterator it = q.begin(); it != q.end(); it++) {
			if (it->touched && (res == q.end() || it->r > res->r)) {
				res = it;
			}
		}
		return res;
	} else {
		return max_element( q.begin(), q.end() );
	}
}

template <typename T> void Queue<T>::simple_add_packet( Packet<T> to_add, bool dosortifneeded ) {
	// preemption for length-aware policies is done later in bulk
	bool wehavespace = lengthaware || (q.size() < B);
	if ( wehavespace ) {
		this->total_admitted++;
		q.push_back(to_add);
		if (dosortifneeded && sorting) doSort();
	} else {
		// preemption for length-aware policies is done in bulk
		if (preempting && !lengthaware) {
			typename vector<Packet<T> >::iterator it = this->get_max_to_preempt();
			if ( it != q.end() && *it > to_add.r * this->beta ) {
				this->total_admitted++;
				q.erase( it );
				q.push_back(to_add);
			}
		}
		if (dosortifneeded && sorting) doSort();
	}
}

template <typename T> void Queue<T>::add_packet( IntPacket i, bool dosortifneeded ) {
	if (multiqueue) {
		mq[i.r - 1].add_packet( i, dosortifneeded );
		return;
	}
	int num_subpackets = 1;
	if (type == "FOPTUW") num_subpackets = i.r;
	if (type == "FOPTUL") num_subpackets = i.l;
	for (int spcount = 0; spcount < num_subpackets; ++spcount) {
		Packet<T> to_add = i;
		if (type == "FOPTUW") { to_add.r = 1; to_add.l =  (i.l / i.r); }
		if (type == "FOPTUL") { to_add.r = (i.r / i.l); to_add.l =  1; }
		this->simple_add_packet(to_add, dosortifneeded);
	}
}

template <typename T> double Queue<T>::totallength() const {
	double res = 0;
	for ( typename vector<Packet<T> >::const_iterator it = q.begin(); it != q.end(); ++it ) {
		res += it->l;
	}
	return res;
}

template class Queue<int>;
template class Queue<float>;


void SharedMultiQueue::print_queue(string pref) {
	ostringstream s; s << pref << "\t[" << type << "]\t";
	for (size_t i=0; i<sq.size(); ++i) {
		s << OUT_SSQ(sq[i]);
	}
	D(s.str());
}

void SharedMultiQueue::add_packets( const vector<IntPacket> & v ) {
	for (size_t i=0; i<v.size(); ++i) {
		add_packet(v[i].r);
	}
}

void SharedMultiQueue::add_packet(int req) {
	int admcheck = check_admission(req);
	if (admcheck >= 0) {
		D("[" << type << "]\tadding packet of req " << req);
		sq[req-1].add_packet();
		++total_admitted;
		if (admcheck >= 1) {
			D("[" << type << "]\tpushing out packet of req " << admcheck << " for a packet of req " << req);
			sq[admcheck-1].drop_packet();
		}
	} else {
		D("[" << type << "]\tdropping packet of req " << req);
	}
}

size_t SharedMultiQueue::size() const {
	size_t res = 0;
	for (int i=0; i<k; ++i) res += sq[i].size;
	return res;
}

void SharedMultiQueue::process(int tick_num) {
	D("[" << type << "]\ttick processing");
	for (size_t j=0; j<C; ++j) {
		for (int i=0; i<k; ++i) {
			if (sq[i].one_tick()) {
				D("[" << type << "]\t\tqueue " << (i+1) << " processed a packet!");
				++total_processed;
			}
		}
	}
}

double harmonic(uint k) {
	double harm = 0;
	for (uint i=1; i<=k; ++i) {
		harm += 1 / (double)i;
	}
	return harm;
}

int NHDTSharedMultiQueue::check_admission(int req) const {
	int cur_size = sq[req-1].size;
	int total = 0; int num = 0;
	for (int i=0; i<k; ++i) {
		if (sq[i].size >= cur_size) {
			total += sq[i].size;
			++num;
		}
	}
	// cout << "[NHDT]\tchecking admission:\tnum=" << num << "\ttotal=" << total << "\tthresh=" << threshold_factor << "\tharm=" << harmonic(num) << endl;
	return ((total < threshold_factor * harmonic(num)) ? 0 : -1 );
};

int LQDSharedMultiQueue::check_admission(int req) const {
	int max_queue = k-1;
	int max_size = sq[k-1].size;
	int total = sq[k-1].size;
	for (int i=k-2; i>=0; --i) {
		total += sq[i].size;
		if (sq[i].size > max_size) {
			max_size = sq[i].size;
			max_queue = i;
		}
	}
	if (total < B) return 0;
	if (max_queue > req) return (max_queue+1);
	return -1;
}

int BPDSharedMultiQueue::check_admission(int req) const {
	int max_nonempty_queue = k-1;
	while (sq[max_nonempty_queue].size == 0 && max_nonempty_queue >= 0) {
		--max_nonempty_queue;
	}
	if (max_nonempty_queue == -1) return 0;
	int total = this->size();
	if (total < B) return 0;
	if (max_nonempty_queue > req) return (max_nonempty_queue+1);
	return -1;
}

int BPD2SharedMultiQueue::check_admission(int req) const {
	int max_nonempty_queue = k-1;
	while (sq[max_nonempty_queue].size < 2 && max_nonempty_queue >= 0) {
		--max_nonempty_queue;
	}
	if (max_nonempty_queue == -1) return 0;
	int total = this->size();
	if (total < B) return 0;
	if (max_nonempty_queue > req) return (max_nonempty_queue+1);
	return -1;
}

int LWDSharedMultiQueue::check_admission(int req) const {
	int total = this->size();
	if (total < B) return 0;
	int max_queue = k-1;
	int max_work = sq[k-1].work();
	for (int i=k-2; i>=0; --i) {
		int cur_work = sq[i].work();
		if (cur_work > max_work) {
			max_work = cur_work;
			max_queue = i;
		}
	}
	if (max_queue > req) return (max_queue+1);
	return -1;
}

