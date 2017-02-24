#pragma once

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/program_options.hpp>

#include "io.hpp"
#include "logging.hpp"
#include "stringutil.hpp"

#include "rules.hpp"


using namespace std;


void do_remove(vector<NSDIRule> & br, const vector<bool> & to_remove) {
	std::vector<NSDIRule> newbr;
	for (uint i=0; i < br.size(); ++i) {
		if (!to_remove[i]) {
			newbr.push_back(br[i]);
		}
	}
	br.swap(newbr);
}

// check for intersections with different actions between low and high with rule r
bool no_intersection_in_the_middle(const vector<NSDIRule> & br, const NSDIRule & r, uint low, uint high) {
	for (uint k=low+1; k<high; ++k) {
		if (r.a != br[k].a && r.intersects(br[k])) {
			return false;
		}
	}
	return true;
}

uint do_forward_subsumption(vector<NSDIRule> & br, uint previous_br_size) {
	vector<bool> to_remove(br.size(), false);
	uint num_subsumed = 0;
	for (uint i=0; i < br.size(); ++i) {
		if (to_remove[i]) continue;
		const NSDIRule & r = br[i];
		// for (uint j=max(previous_br_size, i+1); j < br.size(); ++j) {
		for (uint j=i+1; j < br.size(); ++j) {
			if (to_remove[j]) continue;
			if (r.subsumes(br[j]) == 1) {
				LOG("\t\t" << r.print() << "\t" << br[j].print());
				to_remove[j] = true;
				num_subsumed++;
			}
		}
	}
	do_remove(br, to_remove);
	return num_subsumed;
}

uint do_backward_subsumption(vector<NSDIRule> & br, uint previous_br_size) {
	vector<bool> to_remove(br.size(), false);
	uint num_subsumed = 0;
	for (uint i=0; i < br.size(); ++i) {
		if (to_remove[i]) continue;
		const NSDIRule & r = br[i];
		// for (uint j=max(previous_br_size, i+1); j < br.size(); ++j) {
		for (uint j=i+1; j < br.size(); ++j) {
			// backward subsumption only works if actions are the same
			if (r.a == br[j].a && r.subsumes(br[j]) == -1 && no_intersection_in_the_middle(br, r, i, j)) {
				to_remove[i] = true;
				num_subsumed++;
				break;
			}
		}
	}
	do_remove(br, to_remove);
	return num_subsumed;
}

uint do_resolution(vector<NSDIRule> & br, uint & previous_br_size) {
	vector<bool> to_remove(br.size(), false);
	boost::unordered_map<uint, NSDIRule> to_change;
	int to_be_resolved = -1;
	uint num_resolved = 0;
	for (uint i=0; i < br.size(); ++i) {
		const NSDIRule & r = br[i];
		// for (uint j=max(previous_br_size, i+1); j < br.size(); ++j) {
		for (uint j=i+1; j < br.size(); ++j) {
			if (r.a != br[j].a) continue;
			to_be_resolved = r.can_resolve_by(br[j]);
			if (to_be_resolved >= 0 && no_intersection_in_the_middle(br, r, i, j)) {
				// LOG("\t\t" << i << " " << stringrules[i] << " " << r.print() << "\t" << j << " " << stringrules[j] << " " << br[j].print());
				to_change.insert(make_pair(i, NSDIRule(r, to_be_resolved)));
				to_remove[j] = true;
				num_resolved++;
			}
		}
	}

	// change the rules
	for (auto p : to_change) {
		br[p.first] = p.second;
	}

	do_remove(br, to_remove);
	// here we set which rules we have already processed before
	previous_br_size = br.size();
	return num_resolved;
}

// return a bit which is safe to remove from the relaxed equivalence point of view
vector<bool> find_relaxed_bits(const vector<NSDIRule> & br, uint mode = 0) {
	vector<bool> res(NSDI_BOOL_SIZE, true);
	for (uint i=0; i<br.size(); ++i) {
		for (uint j=i+1; j<br.size(); ++j) {
			br[i].mark_relaxed_intersections(br[j], res, mode);
		}
		if (accumulate(res.begin(), res.end(), 0) == 0) return res;
	}
	return res;
}

vector< vector<bool> > find_relaxed_blocking_subsets(const vector<NSDIRule> & br, const vector<bool> & removed_bits, uint mode = 0) {
	vector< vector<bool> > res = vector< vector<bool> >(br.size(), vector<bool>(NSDI_BOOL_SIZE, false));
	uint min_active_bit = 0;
	for (; min_active_bit < NSDI_BOOL_SIZE; min_active_bit++) {
		if (!removed_bits[min_active_bit]) break;
	}
	uint max_active_bit = NSDI_BOOL_SIZE;
	for (; max_active_bit > 1; max_active_bit--) {
		if (!removed_bits[max_active_bit]) break;
	}
	for (uint i=0; i<br.size(); ++i) {
		for (uint j=i+1; j<br.size(); ++j) {
			br[i].mark_relaxed_intersections_masked(br[j], res[i], removed_bits, min_active_bit, max_active_bit, mode);
		}
	}
	return res;
}

void read_nsdi_file(const string & fname, vector<NSDIRule> & br, boost::unordered_map< string, uint > & actions, vector<string> & action_strings, bool read_binary = false) {
	LOG("Reading file " << fname);
	vector<string> tok;
	vector<uint> max_set;
	process_file_buffered(fname, [&] (uint i, const string & line) {
		boost::algorithm::split(tok, line, boost::is_any_of(" \t"), boost::token_compress_on);
		if (tok.size() < 2) return;
		if (actions.find(tok[1]) == actions.end()) {
			actions[tok[1]] = action_strings.size();
			action_strings.push_back(tok[1]);
		}
		NSDIRule r(tok[0], actions[tok[1]], read_binary);
		// stringrules.push_back(tok[0]);
		br.push_back(r);
	});
	LOG("\tread " << br.size() << " rules...");
}

void apply_binary_rules(vector<NSDIRule> & br) {
	bool changed = true;
	uint iter_count = 0;
	uint previous_br_size = 0;
	while (changed) {
		LOG("Boolean iteration " << iter_count << ":");
		++iter_count;
		changed = false;

		uint num_forward_subsumption = do_forward_subsumption(br, previous_br_size);
		if (num_forward_subsumption > 0) changed = true;
		LOG("\t...forward subsumption removed " << num_forward_subsumption << " rules...");

		uint num_backward_subsumption = do_backward_subsumption(br, previous_br_size);
		if (num_backward_subsumption > 0) changed = true;
		LOG("\t...backward subsumption removed " << num_forward_subsumption << " rules...");

		uint num_resolution = do_resolution(br, previous_br_size);
		if (num_resolution > 0) changed = true;
		LOG("\t...resolution applied " << num_resolution << " times...");

		LOG("\t..." << br.size() << " rules left after the iteration.");
	}
}

uint process_one_optgroup(vector<NSDIRule> & br, std::vector<bool> & subset_d, int num_bits_toremove = -1, uint mode = 0, bool do_binary = true, bool first_iteration = false) {
	if (br.size() == 0) return max(0, num_bits_toremove);
	if (do_binary) {
		apply_binary_rules(br);
	}
	uint min_size = br.size() * NSDI_BOOL_SIZE;
	uint opt_numbits = 0;
	vector<bool> opt_subset_d(br.size(), false);
	std::vector<bool> removed_bits(NSDI_BOOL_SIZE, false);
	bool prev_size_increase = false;
	uint cur_size = br.size() * NSDI_BOOL_SIZE;
	for (int iIter=0; iIter<NSDI_BOOL_SIZE; ++iIter) {
		uint cur_bits_removed = accumulate(removed_bits.begin(), removed_bits.end(), 0);
		uint cur_d_size = accumulate(subset_d.begin(), subset_d.end(), 0);
		int num_bits = cur_bits_removed;

		uint cur_size_new = (NSDI_BOOL_SIZE * br.size() - (br.size() - cur_d_size) * cur_bits_removed );
		LOG( cur_bits_removed << " bits removed. |D| = " << cur_d_size << ", for resulting size " << cur_size_new);
		if (cur_size_new > cur_size && num_bits >= num_bits_toremove ) {
			if (prev_size_increase) {
				LOG("two decreases in size is enough to break");
				break;
			}
			else prev_size_increase = true;
		}
		cur_size = cur_size_new;
		
		if ( ( ((int)num_bits >= num_bits_toremove) && (cur_size < min_size) ) || num_bits == num_bits_toremove ) {
			opt_numbits = num_bits;
			min_size = cur_size;
			for (uint j=0; j<br.size(); ++j) {
				opt_subset_d[j] = subset_d[j];
			}
		}
		vector< vector<bool> > blockers = find_relaxed_blocking_subsets(br, removed_bits, mode);
		vector< uint > zero_blockers;
		uint min_blocker_size = br.size();
		uint min_blocker_ind = -1;
		for (uint i=0; i < NSDI_BOOL_SIZE; ++i) {
			if (removed_bits[i]) continue;
			uint cur_blocker_size = 0;
			for (uint j=0; j<br.size(); ++j) {
				cur_blocker_size += (int)(blockers[j][i]);
			}
			if (cur_blocker_size == 0) {
				zero_blockers.push_back(i);
			}
			if (cur_blocker_size < min_blocker_size) {
				min_blocker_size = cur_blocker_size;
				min_blocker_ind = i;
			}
		}
		if (first_iteration && zero_blockers.size() > 10 && NSDI_BOOL_SIZE - cur_bits_removed > 132 ) {
			uint zero_blockers_to_remove = min( uint(zero_blockers.size() - 5), uint(NSDI_BOOL_SIZE - 32) );
			LOG("Removing " << zero_blockers_to_remove << " bits at once because none of them blocks anything.");
			for (uint i=0; i < zero_blockers_to_remove; ++i) {
				removed_bits[zero_blockers[i]] = true;
			}
		} else {
			// LOG("Next bit to remove is " << min_blocker_ind << " with |D| = " << min_blocker_size);
			removed_bits[min_blocker_ind] = true;
			for (uint i=0; i<br.size(); ++i) {
				if (blockers[i][min_blocker_ind]) {
					subset_d[i] = true;
				}
			}
		}
	}
	subset_d.swap(opt_subset_d);
	return opt_numbits;
}

void process_one_file_optgroup(const string & fname, int num_bits_toremove = -1, uint mode = 0, bool do_first_boolean = true, bool randomize_and_output = false, bool read_binary = false) {
	vector<NSDIRule> br;
	boost::unordered_map< string, uint > actions;
	vector<string> action_strings;
	read_nsdi_file(fname, br, actions, action_strings, read_binary);

	cout << fname << " & " << action_strings.size() << " & " << br.size() << " & " << ((br.size() * NSDI_BOOL_SIZE) / 100) / 10.0;
	if (do_first_boolean) {
		apply_binary_rules(br);
	}
	if (randomize_and_output) {
		LOG("Writing to " + fname + ".randomized.txt...");
		write_to_file(fname + ".randomized.txt", [&] (ofstream & ofs) {
			for (auto r : br) {
				ofs << r.print() << "\t" << r.a << "\n";
			}
		});
		LOG("All done!");
		exit(0);
	}
	double orig_size_kb = ((double)(br.size() * NSDI_BOOL_SIZE) / (double)1000);
	cout << " & " << br.size() << " & " << (int(orig_size_kb * 10) / 10.0);
	uint cur_size = br.size() * NSDI_BOOL_SIZE;
	uint prev_subset_d_numrules = br.size();
	std::vector<bool> subset_d(br.size(), false);
	uint num_removed_bits = process_one_optgroup(br, subset_d, num_bits_toremove, mode, false, true);
	uint cur_subset_d_numrules = accumulate(subset_d.begin(), subset_d.end(), 0);
	uint cur_subset_i_numrules = br.size() - cur_subset_d_numrules;
	cur_size = cur_size - (prev_subset_d_numrules - cur_subset_d_numrules)*num_removed_bits;
	prev_subset_d_numrules = cur_subset_d_numrules;
	LOG("[RESULT] width " << (NSDI_BOOL_SIZE-num_removed_bits) << " bits, " << cur_subset_d_numrules << " rules left, size " << cur_size << ".");
	cout << " & " << (NSDI_BOOL_SIZE-num_removed_bits) << " & " << cur_subset_i_numrules << " & " << cur_subset_d_numrules <<  " & " << (cur_size  / 100) / 10.0;
	uint numGroups = 0;
	for (uint numIter=0; numIter<1000; ++numIter) {
		vector<NSDIRule> cur_br;
		for (uint i=0; i<br.size(); ++i) {
			if (subset_d[i]) cur_br.push_back(br[i]);
		}
		subset_d = std::vector<bool>(cur_br.size(), false);
		br.swap(cur_br);
		apply_binary_rules(br);

		// find the next optgroup
		uint cur_removed_bits = process_one_optgroup(br, subset_d, num_bits_toremove, mode, false, false);

		// update D and I subset sizes
		cur_subset_d_numrules = accumulate(subset_d.begin(), subset_d.end(), 0);
		cur_subset_i_numrules = br.size() - cur_subset_d_numrules;
		cur_size = cur_size - (prev_subset_d_numrules - cur_subset_d_numrules)*num_removed_bits;
		prev_subset_d_numrules = cur_subset_d_numrules;

		LOG("[RESULT] Group " << (numIter+2) << " done, width " << (NSDI_BOOL_SIZE-cur_removed_bits) << " bits, current |I|=" << cur_subset_i_numrules << ", |D|=" << cur_subset_d_numrules << ", total size=" << (cur_size  / 100) / 10.0);
		if (numIter < 3) cout << " & " << (NSDI_BOOL_SIZE-cur_removed_bits) << " & " << cur_subset_i_numrules << " & " << cur_subset_d_numrules << " & " << (cur_size  / 100) / 10.0;
		if (cur_subset_d_numrules == 0) {
			if (numGroups == 0) numGroups = numIter+1;
			if (numIter >= 3) break;
		}
	}
	double cur_size_kb = (double)cur_size / (double)1000.0;
	cout << " & {\\bf " << (numGroups+1) << "} & {\\bf " << (int(cur_size_kb * 10) / 10.0) << "} & {\\bf " << (int)(1000.0 * cur_size_kb / orig_size_kb) / 1000.0 << "}\\\\\n";
}

