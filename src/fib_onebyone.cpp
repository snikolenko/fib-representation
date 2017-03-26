#include "reduce_fib.hpp"

using namespace std;

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
	std::vector<int> group_ids;
} AllGroupsResult;

void print_AGR(const AllGroupsResult & res) {
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
		cout << "Processing " << cur_br.size() << " rules, removing " << bits_toremove << " bits" << "...\n";
		OneOptGroupResult ogres = process_one_optgroup_L(cur_br, bits_toremove, mode, false, true);
		cout << "Group " << iGroup << ":\n";
		print_OGR(ogres, cur_br);

		// add the current group to AllGroupsResult
		res.group_bits.push_back( ogres.group_bits );
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

uint add_rule_to_groups(const NSDIRule & r, vector<NSDIRule> & br, vector<NSDIRule> & D_set, AllGroupsResult & agr, uint bits_toremove, uint mode = 1) {
	// check if we can add the new rule to one of the groups
	// cout << "Adding rule\n" << r.print() << endl;
	for (uint iGroup=0; iGroup < agr.group_bits.size(); ++iGroup) {
		// cout << "Against group " << iGroup << ":\n" << r.print(agr.group_bits[iGroup]) << endl;
		// check if we can add to this group
		bool can_add = true;
		for (uint i=0; i<agr.group_ids.size(); ++i) {
			// check intersections with every rule in the group w.r.t. the current group bits
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
			return 0;
		}
	}
	// if couldn't add, try to leave in D
	if (D_set.size() < D_size) {
		D_set.push_back(r);
		return 0;
	}
	// if couldn't, recompute everything
	cout << "Recomputing everything! D_set size = " << D_set.size() << endl;
	br.push_back(r);
	br.insert(br.end(), D_set.begin(), D_set.end());
	agr = recompute_groups(br, bits_toremove);
	print_AGR(agr);
	return 1;
}

int main(int argc, char* argv[]) {
	srand(42);

	po::options_description desc("Program options");
	desc.add_options()
	    ("help,?", "produce help message")
	    ("max_width,w", po::value<int>(&max_width)->default_value(-1), "max group width in bits")
	    ("mode,m", po::value<int>(&mode)->default_value(1), "mode")
	    ("d_size,d", po::value<uint>(&D_size)->default_value(0), "size of D")
	    ("input,i", po::value<string>(&in_fname)->default_value(""), "input file")
	    ("ingress,r", po::value<string>(&ingress_fname)->default_value(""), "input file for the ingress part of the classifier")
	    ("egress0,0", po::value<string>(&egress0_fname)->default_value(""), "input file for the egress part 0 of the classifier")
	    ("egress1,1", po::value<string>(&egress1_fname)->default_value(""), "input file for the egress part 1 of the classifier")
	    ("nobinary,b", "do not apply binary reduction before first iteration")
	    ("addrandom,a", "add random bits and output the result")
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

	int bits_toremove = -1;
	if (max_width > -1) {
		bits_toremove = NSDI_BOOL_SIZE - max_width;
	}

	vector<NSDIRule> input_br;
	boost::unordered_map< string, uint > actions;
	vector<string> action_strings;
	read_nsdi_file(in_fname, input_br, actions, action_strings);

	vector<NSDIRule> br;
	for (uint i=input_br.size() - 100; i<input_br.size(); ++i) {
		br.push_back(input_br[i]);
	}

	AllGroupsResult agr = recompute_groups(br, bits_toremove);
	print_AGR(agr);

	uint num_recomputations = 0;
	vector<NSDIRule> D_set;
	for (int i=input_br.size() - 101; i >= 0; --i) {
		cout << i << endl;
		uint res_add = add_rule_to_groups(input_br[i], br, D_set, agr, bits_toremove, mode);
		num_recomputations = num_recomputations + res_add;
		cout << "\t\tcurrent recomputations: " << num_recomputations << endl;
	}
	cout << "Total recomputations: " << num_recomputations << endl;
	exit(0);
}

