/*
 * queue.hpp
 *
 *  Created on: Jun 8, 2012
 *      Author: snikolenko
 */

#ifndef QUEUE_HPP_
#define QUEUE_HPP_

#include <string>
#include <vector>
#include <iostream>
 #include <math.h>

#include "packet.hpp"

// #define DEBUG
#ifdef DEBUG
	#define D(a) cout << a << endl
#else
	#define D(a)
#endif

#define EPSILON 0.001

using namespace std;
typedef unsigned int uint;

// forward declaration
template <typename T> struct t_compare_queues_by_length;

template <typename T> class Queue {
public:
	Queue(string iType, uint ik, uint iB, uint iC, bool iMultiqueue) : type(iType), B(iB), C(iC), L(1), k(ik), beta(1),
			total_processed(0), total_processed_length(0), total_admitted(0), total_delay(0), total_squared_delay(0),
			multiqueue(iMultiqueue), recycling(false), preempting(false), sorting(false), lengthaware(false),
			fractional(false), lazy(false), reversing(false), leavepartiallyprocessed(false), sendingout(false) {
		init();
	}
	Queue(string iType, uint ik, uint iB, uint iC, bool iMultiqueue, double bt) : type(iType), B(iB), C(iC), L(1), k(ik), beta(bt),
			total_processed(0), total_processed_length(0), total_admitted(0), total_delay(0), total_squared_delay(0),
			multiqueue(iMultiqueue), recycling(false), preempting(false), sorting(false), lengthaware(false),
			fractional(false), lazy(false), reversing(false), leavepartiallyprocessed(false), sendingout(false) {
		init();
	}
	Queue(string iType, uint ik, uint iB, uint iC, uint iL, bool iMultiqueue) : type(iType), B(iB), C(iC), L(iL), k(ik), beta(1),
			total_processed(0), total_processed_length(0), total_admitted(0), total_delay(0), total_squared_delay(0),
			multiqueue(iMultiqueue), recycling(false), preempting(false), sorting(false), lengthaware(false),
			fractional(false), lazy(false), reversing(false), leavepartiallyprocessed(false), sendingout(false) {
		init();
	}

	string type;
	size_t B;
	size_t C;
	size_t L;
	size_t k;

	double beta;

	T total_processed;
	double total_processed_length;
	int total_admitted;
	unsigned int total_delay;
	unsigned int total_squared_delay;

	bool multiqueue;
	bool recycling;
	bool preempting;
	bool sorting;
	bool lengthaware;
	bool fractional;
	bool lazy;
	bool reversing;
	bool leavepartiallyprocessed;

	bool sendingout;

	vector< Packet<T> > q;
	vector< Queue > mq;
	size_t curq;

	t_compare_queues_by_length<T> compare_queues_by_length;

	bool has_packets() const { return !q.empty(); }
	const Packet<T> hol() const { return (q.empty() ? Packet<T>(-1, -1) : q[0]); }
	typename vector<Packet<T> >::iterator get_max_to_preempt();
	virtual size_t size() const { return q.size(); }
	double totallength() const;
	void increment_counters(int length, int delay, int squared_delay);

	void print_queue(string pref);

	void process_head_packet(int tick_num);
	void remove_head_packet(int tick_num);
	void simple_add_packet( Packet<T> to_add, bool dosortifneeded );
	void add_packet( IntPacket i, bool dosortifneeded = true );
	void add_packets( const vector<IntPacket> & v );

	virtual void process(int tick_num);
	void process_multiqueue(int tick_num);
	void process_subqueue(int tick_num, int i);
	void doSort();
	void init();

	int get_next_queue(size_t cur);

	void recycle();
};

#define OUT_SSQ(a) "[" << (a).k << ":" << (a).size << ":" << (a).hol << "] "
#define OUT_THIS_SSQ "[" << k << ":" << size << ":" << hol << "] "
class SharedSubQueue {
public:
	size_t size;
	size_t hol;
	int k;
	int work() const {
		return (size == 0 ? 0 : ((size-1) * k + hol));
	}
	SharedSubQueue(int ik) : size(0), hol(0), k(ik) {}
	void drop_packet() {
		if (size == 0) {
			cout << "\t\t\tERROR: dropping packet from empty queue" << endl;
		} else {
			size--;
			if (size == 0) hol = 0;
		}
	}
	void add_packet() {
		if (size == 0) hol = k;
		size++;
	}
	bool one_tick() {
		D("\t\t\tone tick for " << OUT_THIS_SSQ);
		if (hol > 1) {
			hol--;
		} else if (hol == 1) {
			if (size > 0) {
				size--;
				if (size > 0) hol = k;
				else hol = 0;
				return true;
			}
		}
		return false;
	}
};

class SharedMultiQueue {
public:
	string type;
	size_t k, B, C, L;
	double beta;
	int total_processed, total_admitted, total_processed_length, total_delay, total_squared_delay;

	// vector of subqueues; a subqueue is its size + remaining 
	vector< SharedSubQueue > sq;
	SharedMultiQueue(string iType, uint ik, uint iB, uint iC): type(iType), k(ik), B(iB), C(iC), L(1), beta(1.0),
			total_processed(0), total_admitted(0), total_processed_length(0), total_delay(0), total_squared_delay(0) {
		sq.clear();
		for (size_t i=1; i<=ik; ++i) {
			sq.push_back(i);
		}
		for (size_t i=0; i<ik; ++i) {
			D(OUT_SSQ(sq[i]));
		}
	}
	virtual int check_admission(int req) const = 0;
	void add_packets( const vector<IntPacket> & v );
	virtual void add_packet(int req);
	virtual size_t size() const;
	virtual void process(int tick_num);
	void print_queue(string pref);
};

double harmonic(uint k);

class NHSTSharedMultiQueue : public SharedMultiQueue {
public:
	vector<double> threshold;
	NHSTSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("NHST", k, B, C) {
		double numerator = (double)B / harmonic(k);
		for (size_t i=1; i<=k; ++i) {
			threshold.push_back(numerator / i);
		}
	}
	virtual int check_admission(int req) const {return (sq[req-1].size < threshold[req-1]) ? 0 : -1; };
};

class NESTSharedMultiQueue : public SharedMultiQueue {
public:
	vector<double> threshold;
	NESTSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("NEST", k, B, C) {
		for (size_t i=1; i<=k; ++i) {
			threshold.push_back(B / k);
		}
	}
	virtual int check_admission(int req) const {return (sq[req-1].size < threshold[req-1]) ? 0 : -1; };
};

class NHDTSharedMultiQueue : public SharedMultiQueue {
public:
	double threshold_factor;
	NHDTSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("NHDT", k, B, C) {
		threshold_factor = B / harmonic(k);
	}
	virtual int check_admission(int req) const;
};

class LQDSharedMultiQueue : public SharedMultiQueue {
public:
	LQDSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("LQD", k, B, C) { }
	virtual int check_admission(int req) const;
};

class BPDSharedMultiQueue : public SharedMultiQueue {
public:
	BPDSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("BPD", k, B, C) { }
	virtual int check_admission(int req) const;
};

class BPD2SharedMultiQueue : public SharedMultiQueue {
public:
	BPD2SharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("BPD2", k, B, C) { }
	virtual int check_admission(int req) const;
};

class LWDSharedMultiQueue : public SharedMultiQueue {
public:
	LWDSharedMultiQueue(uint k, uint B, uint C) : SharedMultiQueue("LWD", k, B, C) { }
	virtual int check_admission(int req) const;
};

template <typename T> struct t_compare_queues_by_length {
	bool operator () ( const Queue<T> & q1, const Queue<T> & q2 ) {
		return q1.size() < q2.size();
	}
};

class QueueContainer {
public:
	virtual const string & type() const = 0;
	virtual size_t B() const = 0;
	virtual size_t C() const = 0;
	virtual size_t L() const = 0;
	virtual size_t k() const = 0;
	virtual double beta() const = 0;

	virtual unsigned int total_delay() const = 0;
	virtual unsigned int total_squared_delay() const = 0;
	virtual unsigned int total_processed() const = 0;
	virtual unsigned int total_processed_length() const = 0;
	virtual unsigned int final_total_processed() const = 0;
	virtual unsigned int total_admitted() const = 0;
	virtual double final_total_processed_length() const = 0;

	virtual void add_packets( const vector<IntPacket> & v ) = 0;
	virtual void process(int tick_num) = 0;
	virtual void print_queue(string pref) = 0;
};

class IntQueueContainer : public QueueContainer {
public:
	IntQueueContainer(string iType, uint ik, uint iB, uint iC, bool iMultiqueue) : q(iType, ik, iB, iC, iMultiqueue) {};
	IntQueueContainer(string iType, uint ik, uint iB, uint iC, bool iMultiqueue, double bt) : q(iType, ik, iB, iC, iMultiqueue, bt) {};
	IntQueueContainer(string iType, uint ik, uint iB, uint iC, uint iL, bool iMultiqueue) : q(iType, ik, iB, iC, iL, iMultiqueue) {};

	virtual const string & type() const { return q.type; }
	virtual size_t B() const { return q.B; }
	virtual size_t C() const { return q.C; }
	virtual size_t L() const { return q.L; }
	virtual size_t k() const { return q.k; }
	virtual double beta() const { return q.beta; }

	virtual unsigned int total_delay() const { return q.total_delay; }
	virtual unsigned int total_squared_delay() const { return q.total_squared_delay; }
	virtual unsigned int total_processed() const { return q.total_processed; }
	virtual unsigned int total_admitted() const { return q.total_admitted; }
	virtual unsigned int total_processed_length() const { return q.total_processed_length; }
	virtual unsigned int final_total_processed() const { return q.total_processed + q.size(); }
	virtual double final_total_processed_length() const { return q.total_processed_length + q.totallength(); }

	virtual void add_packets( const vector<IntPacket> & v ) { q.add_packets(v); }
	virtual void process(int tick_num) { q.process(tick_num); }
	virtual void print_queue(string pref) { q.print_queue(pref); }
private:
	Queue<int> q;
};

class FloatQueueContainer : public QueueContainer {
public:
	FloatQueueContainer(string iType, uint ik, uint iB, uint iC, bool iMultiqueue) : q(iType, ik, iB, iC, iMultiqueue) {};
	FloatQueueContainer(string iType, uint ik, uint iB, uint iC, bool iMultiqueue, double bt) : q(iType, ik, iB, iC, iMultiqueue, bt) {};
	FloatQueueContainer(string iType, uint ik, uint iB, uint iC, uint iL, bool iMultiqueue) : q(iType, ik, iB, iC, iL, iMultiqueue) {};

	virtual const string & type() const { return q.type; }
	virtual size_t B() const { return q.B; }
	virtual size_t C() const { return q.C; }
	virtual size_t L() const { return q.L; }
	virtual size_t k() const { return q.k; }
	virtual double beta() const { return q.beta; }

	virtual unsigned int total_delay() const { return q.total_delay; }
	virtual unsigned int total_squared_delay() const { return q.total_squared_delay; }
	virtual unsigned int total_processed() const { return q.total_processed; }
	virtual unsigned int total_admitted() const { return q.total_admitted; }
	virtual unsigned int total_processed_length() const { return q.total_processed_length; }
	virtual unsigned int final_total_processed() const { return q.total_processed + q.size(); }
	virtual double final_total_processed_length() const { return q.total_processed_length + q.totallength(); }

	virtual void add_packets( const vector<IntPacket> & v ) { q.add_packets(v); }
	virtual void process(int tick_num) { q.process(tick_num); }
	virtual void print_queue(string pref) { q.print_queue(pref); }
private:
	Queue<float> q;
};

class SMQContainer : public QueueContainer {
public:
	SMQContainer(SharedMultiQueue *qPtr) : q(qPtr) {};
	
	virtual const string & type() const { return q->type; }
	virtual size_t B() const { return q->B; }
	virtual size_t C() const { return q->C; }
	virtual size_t L() const { return q->L; }
	virtual size_t k() const { return q->k; }
	virtual double beta() const { return q->beta; }

	virtual unsigned int total_delay() const { return q->total_delay; }
	virtual unsigned int total_squared_delay() const { return q->total_squared_delay; }
	virtual unsigned int total_processed() const { return q->total_processed; }
	virtual unsigned int total_admitted() const { return q->total_admitted; }
	virtual unsigned int total_processed_length() const { return q->total_processed_length; }
	virtual unsigned int final_total_processed() const { return q->total_processed + q->size(); }
	virtual double final_total_processed_length() const { return q->total_processed_length + q->size(); }

	virtual void add_packets( const vector<IntPacket> & v ) { q->add_packets(v); }
	virtual void process(int tick_num) { q->process(tick_num); }
	virtual void print_queue(string pref) { q->print_queue(pref); }
private:
	SharedMultiQueue *q;
};

#endif /* QUEUE_HPP_ */
