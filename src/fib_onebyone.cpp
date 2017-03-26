#include "reduce_fib.hpp"

using namespace std;

#include <algorithm>

namespace po = boost::program_options;

int max_width = -1;

string in_fname = "";
string ingress_fname = "";
string egress0_fname = "";
string egress1_fname = "";

bool do_first_boolean = true;
bool randomize_and_output = false;
bool read_binary = false;

uint D_size = 0;
uint beta = 0;
uint random_size = 50000;
int mode = 1;

void print_result(const string & key, const string & value) {
	cout << key << "\t" << value << "\n";
	cout.flush();
}

void print_result(const string & key, int value) {
	cout << key << "\t" << value << "\n";
	cout.flush();
}

void print_OGR(const OneOptGroupResult & res, const vector<NSDIRule> & br) {
	cout << "bits:\t";
	for (uint i=0; i<NSDI_BOOL_SIZE; ++i) {
		cout << res.group_bits[i];
	}
	cout << endl;

	cout << "D/I:\t";
	for (uint i=0; i<br.size(); ++i) {
		if (res.subset_d[i] && res.subset_i[i]) {
			cout << 'X';
		} else if (res.subset_d[i]) {
			cout << "D";
		} else if (res.subset_i[i]) {
			cout << "_";
		} else {
			cout << "?";
		}
	}
	cout << endl;
}

OneOptGroupResult process_one_optgroup_L(vector<NSDIRule> & br, int num_bits_toremove = -1, uint mode = 0, bool do_binary = true, bool first_iteration = false) {
	OneOptGroupResult res;
	if (br.size() == 0) {
		res.opt_numbits = max(0, num_bits_toremove);
		return res;
	}
	if (do_binary) {
		apply_binary_rules(br);
	}
	uint min_size = br.size() * NSDI_BOOL_SIZE;
	uint opt_numbits = 0;
	res.subset_d = vector<bool>(br.size(), false);
	res.subset_i = vector<bool>(br.size(), true);
	vector<bool> opt_subset_d = vector<bool>(br.size(), false);
	vector<bool> removed_bits(NSDI_BOOL_SIZE, false);
	// bool prev_size_increase = false;
	uint cur_size = br.size() * NSDI_BOOL_SIZE;
	for (int iIter=0; iIter<NSDI_BOOL_SIZE; ++iIter) {
		uint cur_bits_removed = accumulate(removed_bits.begin(), removed_bits.end(), 0);
		uint cur_d_size = accumulate(res.subset_d.begin(), res.subset_d.end(), 0);
		int num_bits = cur_bits_removed;

		cur_size = (NSDI_BOOL_SIZE * br.size() - (br.size() - cur_d_size) * cur_bits_removed );
		
		if ( ( ((int)num_bits >= num_bits_toremove) && (cur_size < min_size) ) || num_bits == num_bits_toremove ) {
			opt_numbits = num_bits;
			min_size = cur_size;
			for (uint j=0; j<br.size(); ++j) {
				opt_subset_d[j] = res.subset_d[j];
			}
			if (cur_bits_removed >= num_bits_toremove) {
				LOG("Removed " << num_bits_toremove << " bits, exiting.");
				break;
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
			LOG("Next bit to remove is " << min_blocker_ind << " with |D| = " << min_blocker_size);
			removed_bits[min_blocker_ind] = true;
			for (uint i=0; i<br.size(); ++i) {
				if (blockers[i][min_blocker_ind]) {
					res.subset_d[i] = true;
				}
			}
		}
	}
	res.subset_d.swap(opt_subset_d);
	for (uint i=0; i<br.size(); ++i) {
		res.subset_i[i] = (res.subset_i[i] && (!res.subset_d[i]));
	}
	res.opt_numbits = opt_numbits;
	res.group_bits = std::vector<bool>(NSDI_BOOL_SIZE, true);
	for (uint i=0; i<NSDI_BOOL_SIZE; ++i) {
		res.group_bits[i] = !removed_bits[i];
	}
	return res;
}



typedef struct  {
	vector<uint> widths, num_rules;
	std::vector< std::vector<bool> > group_bits;
	// std::vector< std::vector<uint> > group_bitindices;
	std::vector<int> group_ids;
} AllGroupsResult;

void print_AGR(const AllGroupsResult & res, bool print_groupids = false) {
	cout << "#rules:\t";
	for (uint i=0; i<res.num_rules.size(); ++i) {
		cout << res.num_rules[i] << " ";
	}
	cout << endl;

	cout << "widths:\t";
	for (uint i=0; i<res.widths.size(); ++i) {
		cout << res.widths[i] << " ";
	}
	cout << endl;

	cout << "bits:\n";
	for (uint iGroup=0; iGroup<res.group_bits.size(); ++iGroup) {
		for (uint i=0; i<NSDI_BOOL_SIZE; ++i) {
			cout << res.group_bits[iGroup][i];
		}
		cout << endl;
	}

	if (print_groupids) {
		cout << "gr_ids:\t";
		for (uint i=0; i<res.group_ids.size(); ++i) {
			cout << res.group_ids[i];
		}
		cout << endl;
	}
}

AllGroupsResult	recompute_groups(const vector<NSDIRule> & input_br, uint bits_toremove, int mode = 1) {
	// init br, agr, rule_indices
	vector<NSDIRule> cur_br(input_br.begin(), input_br.end());
	AllGroupsResult	res;
	res.group_ids = vector<int>(cur_br.size(), -1);
	vector<uint> rule_indices(cur_br.size());
	for (uint i=0; i<cur_br.size(); ++i) {
		rule_indices[i] = i;
	}

	for (uint iGroup=0; iGroup<1000; ++iGroup) {
		// cout << "Processing " << cur_br.size() << " rules, removing " << bits_toremove << " bits" << "...\n";
		OneOptGroupResult ogres = process_one_optgroup_L(cur_br, bits_toremove, mode, false, true);
		// cout << "Group " << iGroup << ":\n";
		// print_OGR(ogres, cur_br);

		// add the current group to AllGroupsResult
		res.group_bits.push_back( ogres.group_bits );
		// vector<uint> cur_bitindices;
		// for (uint i=0; i<NSDI_BOOL_SIZE; ++i) {
		// 	if (ogres.group_bits[i]) {
		// 		cur_bitindices.push_back( i );
		// 	}
		// }
		// res.group_bitindices.push_back( cur_bitindices );
		res.widths.push_back( accumulate(ogres.group_bits.begin(), ogres.group_bits.end(), 0) );
		res.num_rules.push_back( accumulate(ogres.subset_i.begin(), ogres.subset_i.end(), 0) );
		for (uint i=0; i<cur_br.size(); ++i) {
			if (ogres.subset_i[i]) {
				res.group_ids[ rule_indices[i] ] = iGroup;
			}
		}
		// recompute indices for next iteration
		vector<NSDIRule> new_br;
		rule_indices.clear();
		for (uint i=0; i<cur_br.size(); ++i) {
			if (ogres.subset_d[i]) {
				rule_indices.push_back(i);
				new_br.push_back(cur_br[i]);
			}
		}
		cur_br.swap(new_br);
		// exit if no more rules
		if (cur_br.empty()) {
			break;
		}
	}
	return res;
}

// uint 

uint add_rule_to_groups(const NSDIRule & r, vector<NSDIRule> & br, vector<NSDIRule> & D_set, AllGroupsResult & agr, uint bits_toremove, uint mode = 1, uint max_D_size=0) {
	// check if we can add the new rule to one of the groups
	// cout << "Adding rule\n" << r.print() << endl;
	for (uint iGroup=0; iGroup < agr.group_bits.size(); ++iGroup) {
		// cout << "Against group " << iGroup << ":\n" << r.print(agr.group_bits[iGroup]) << endl;
		// check if we can add to this group
		bool can_add = true;
		for (uint i=0; i<agr.group_ids.size(); ++i) {
			// check intersections with every rule in the group w.r.t. the current group bits
			// if ( (agr.group_ids[i] == (int)iGroup) && r.intersects_in_mode_masked_above(br[i], agr.group_bitindices[iGroup], mode) ) {
			if ( (agr.group_ids[i] == (int)iGroup) && r.intersects_in_mode_masked_above(br[i], agr.group_bits[iGroup], mode) ) {
				// cout << "Group " << iGroup << ": intersection with rule\n" << br[i].print(agr.group_bits[iGroup]) << endl;
				can_add = false;
				break;
			}
		}
		if (can_add) {
			// add to this group, update agr, all is well, return false
			br.push_back(r);
			agr.num_rules[iGroup] += 1;
			agr.group_ids.push_back(iGroup);
			// cout << "\t\tadding to group " << iGroup << "...";
			// print_AGR(agr, true);
			return 0;
		}
	}
	// if couldn't add, try to leave in D
	if (D_set.size() < max_D_size) {
		D_set.push_back(r);
		return 0;
	}
	// if couldn't, recompute everything
	// cout << "Recomputing everything! D_set size = " << D_set.size() << endl;
	br.push_back(r);
	br.insert(br.end(), D_set.begin(), D_set.end());
	D_set.clear();
	agr = recompute_groups(br, bits_toremove);
	// print_AGR(agr);
	return 1;
}

uint get_recomputations(const vector<NSDIRule> & input_br, uint beta, uint max_width, uint max_D_size) {
	int bits_toremove = NSDI_BOOL_SIZE - max_width;

	vector<NSDIRule> br;
	for (uint i=input_br.size() - 100; i<input_br.size(); ++i) {
		br.push_back(input_br[i]);
	}

	// cout << "Beta=" << beta << "\tmax.wid=" << max_width << "\tmax|D|=" << max_D_size << endl;

	AllGroupsResult agr = recompute_groups(br, bits_toremove);
	// print_AGR(agr);
	if (agr.group_bits.size() < beta) {
		for (uint iGroup=agr.group_bits.size(); iGroup < beta; ++iGroup) {
			vector<bool> new_groupbits(NSDI_BOOL_SIZE, false);
			vector<uint> new_groupbitindices;
			for (uint i=0; i<max_width; ++i) {
				new_groupbits[i] = true;
				// new_groupbitindices.push_back(i);
			}
			agr.group_bits.push_back(new_groupbits);
			// agr.group_bitindices.push_back(new_groupbitindices);
			agr.widths.push_back(max_width);
			agr.num_rules.push_back(0);
		}
	}
	// print_AGR(agr);

	uint num_recomputations = 0;
	vector<NSDIRule> D_set;
	for (int i=input_br.size() - 101; i >= 0; --i) {
		if (i % 10000 == 0) {
			// LLOG(i << "  |D|=" << D_set.size());
		}
		uint res_add = add_rule_to_groups(input_br[i], br, D_set, agr, bits_toremove, mode, max_D_size);
		num_recomputations = num_recomputations + res_add;
		if (res_add == 1) {
			// cout << "recomp. # " << num_recomputations << " at step " << i << endl;
		}
		// cout << "\t\tcurrent recomputations: " << num_recomputations << endl;
	}
	return num_recomputations;
}

int main(int argc, char* argv[]) {
	srand(42);

	po::options_description desc("Program options");
	desc.add_options()
	    ("help,?", "produce help message")
	    ("max_width,w", po::value<int>(&max_width)->default_value(-1), "max group width in bits")
	    ("mode,m", po::value<int>(&mode)->default_value(1), "mode")
	    ("beta,t", po::value<uint>(&beta)->default_value(1), "beta")
	    ("d_size,d", po::value<uint>(&D_size)->default_value(0), "size of D")
	    ("input,i", po::value<string>(&in_fname)->default_value(""), "input file")
	    ("ingress,r", po::value<string>(&ingress_fname)->default_value(""), "input file for the ingress part of the classifier")
	    ("egress0,0", po::value<string>(&egress0_fname)->default_value(""), "input file for the egress part 0 of the classifier")
	    ("egress1,1", po::value<string>(&egress1_fname)->default_value(""), "input file for the egress part 1 of the classifier")
	    ("nobinary,b", "do not apply binary reduction before first iteration")
	    ("addrandom,a", "add random bits and output the result")
	    ("tablerow,2", "do a whole table row for the experiments")
	    ("tablerow2,3", "do a whole table row for the experiments -- 2")
	    ("tablerow3,4", "do a whole table row for the experiments -- 3")
	    ("random,9", "do randomized trials for 50K each")
	    ("randsize,8", po::value<uint>(&random_size)->default_value(50000), "beta")
	    ("readbinary,c", "read binary result immediately")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	if (in_fname == "") {
	    cout << "Please specify input file name.\n";
	    return 1;
	}

	vector<NSDIRule> input_br;
	boost::unordered_map< string, uint > actions;
	vector<string> action_strings;
	read_nsdi_file(in_fname, input_br, actions, action_strings);

	cout << in_fname;
	flush(cout);

	if (vm.count("random")) {
		for (uint cur_beta : { beta }) {
			for (uint max_width : { 13, 16, 24 }) {
				for (uint D_size : { 500, 1000, 2000 }) {
					vector<uint> v;
					for (uint iRun=0; iRun<6; iRun++) {
						random_shuffle( input_br.begin(), input_br.end() );
						uint num_recomputations;
						if (NSDI_BOOL_SIZE == 128) {
							num_recomputations = get_recomputations(input_br, cur_beta, max_width, D_size);
						} else {
							vector<NSDIRule> current_br( input_br.begin(), input_br.begin() + random_size );
							num_recomputations = get_recomputations(current_br, beta, max_width, D_size);
						}
						v.push_back(num_recomputations);
						// cout << " " << num_recomputations << endl;
					}

					double sum = std::accumulate(v.begin(), v.end(), 0.0);
					double mean = sum / v.size();

					double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
					double stddev = std::sqrt(sq_sum / v.size() - mean * mean);

					cout << " & $" << mean << "\\pm " << setprecision(2) << stddev << "$";
					flush(cout);
				}
			}
		}

		exit(0);
	}

	if (vm.count("tablerow")) {
		for (uint beta : { 0, 2, 4 }) {
			for (uint max_width : { 13, 16, 24 }) {
				for (uint D_size : { 0, 20, 50 }) {
					uint num_recomputations = get_recomputations(input_br, beta, max_width, D_size);
					cout << " & " << num_recomputations;
					flush(cout);
				}
			}
		}
	} else if (vm.count("tablerow2")) {
		for (uint beta : { 4 }) {
			for (uint max_width : { 13, 16, 24 }) {
				for (uint D_size : { 1000, 2000, 5000 }) {
					uint num_recomputations = get_recomputations(input_br, beta, max_width, D_size);
					cout << " & " << num_recomputations;
					flush(cout);
				}
			}
		}
	} else if (vm.count("tablerow3")) {
		for (uint cur_beta : { beta }) {
			for (uint max_width : { 13, 16, 24 }) {
				for (uint D_size : { 1000, 2000, 5000 }) {
					uint num_recomputations = get_recomputations(input_br, cur_beta, max_width, D_size);
					cout << " & " << num_recomputations;
					flush(cout);
				}
			}
		}
	} else if (vm.count("tablerow3")) {
		for (uint beta : { 8 }) {
			for (uint max_width : { 13, 16, 24 }) {
				for (uint D_size : { 1000, 2000, 5000 }) {
					uint num_recomputations = get_recomputations(input_br, beta, max_width, D_size);
					cout << " & " << num_recomputations;
					flush(cout);
				}
			}
		}
	} else {
		uint num_recomputations = get_recomputations(input_br, beta, max_width, D_size);
		cout << num_recomputations;
	}

	cout << "\\\\\n";

	exit(0);
}

