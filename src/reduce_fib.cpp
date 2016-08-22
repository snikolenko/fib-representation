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

void print_result(const string & key, const string & value) {
	cout << key << "\t" << value << "\n";
	cout.flush();
}

void print_result(const string & key, int value) {
	cout << key << "\t" << value << "\n";
	cout.flush();
}

void print_result_float(const string & key, double value) {
	cout << key << "\t" << std::fixed << std::setprecision(1) << value << "\n";
	cout.flush();
}

typedef struct  {
	uint num_rules_final, num_rules_1group, width_1group, I_1group, D_1group, num_groups;
	float size_final, size_1group, egress_fpcheck_size;
} OneModeResult;


void print_OMR(const string & prefix, const OneModeResult & omr) {
	print_result(prefix + " 1gr wid", omr.width_1group);
	print_result(prefix + " 1gr I", omr.I_1group);
	print_result(prefix + " 1gr D", omr.D_1group);
	print_result_float(prefix + " 1gr size", omr.size_1group);
	print_result(prefix + " #gr", omr.num_groups);
	print_result_float(prefix + " size", omr.size_final);
	print_result_float(prefix + " egress size", omr.egress_fpcheck_size);
}

OneModeResult process_one_mode(const vector<NSDIRule> & input_br, uint mode, int num_bits_toremove) {
	vector<NSDIRule> br;
	OneModeResult res = {};
	for (uint i=0; i<input_br.size(); ++i) {
		br.push_back(input_br[i]);
	}

	uint cur_size = br.size() * NSDI_BOOL_SIZE;
	uint prev_subset_d_numrules = br.size();
	std::vector<bool> subset_d(br.size(), false);

	uint egress_fpcheck_size = 0;
	uint num_removed_bits = process_one_optgroup(br, subset_d, num_bits_toremove, mode, false, true);
	uint cur_subset_d_numrules = accumulate(subset_d.begin(), subset_d.end(), 0);
	uint cur_subset_i_numrules = br.size() - cur_subset_d_numrules;
	egress_fpcheck_size += (NSDI_BOOL_SIZE - num_removed_bits) * cur_subset_i_numrules;
	cur_size = cur_size - (prev_subset_d_numrules - cur_subset_d_numrules)*num_removed_bits;
	prev_subset_d_numrules = cur_subset_d_numrules;
	uint numGroups = 0;
	res.width_1group = NSDI_BOOL_SIZE-num_removed_bits;
	res.I_1group = cur_subset_i_numrules;
	res.D_1group = cur_subset_d_numrules;
	res.size_1group = egress_fpcheck_size / 1000.;
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
		egress_fpcheck_size += (NSDI_BOOL_SIZE - cur_removed_bits) * cur_subset_i_numrules;
		cur_size = cur_size - (prev_subset_d_numrules - cur_subset_d_numrules)*num_removed_bits;
		prev_subset_d_numrules = cur_subset_d_numrules;

		if (cur_subset_d_numrules == 0) {
			if (numGroups == 0) numGroups = numIter+1;
			break;
		}
	}
	double cur_size_kb = (double)cur_size / (double)1000.0;
	res.num_groups = numGroups+1;
	res.size_final = cur_size_kb;
	res.egress_fpcheck_size = egress_fpcheck_size / 1000.;
	return res;
}

float get_size_round(uint num_rules, uint rule_width) {
	return (int(((double)(num_rules * rule_width) / (double)1000) * 10) / 10.0);
}

void process_file_allmodes(const string & fname, int num_bits_toremove = -1) {
	vector<NSDIRule> input_br;
	boost::unordered_map< string, uint > actions;
	vector<string> action_strings;
	read_nsdi_file(fname, input_br, actions, action_strings);

	// print original sizes
	float size_original = get_size_round(input_br.size(), NSDI_BOOL_SIZE);
	print_result("Original actions", action_strings.size());
	print_result("Original rules", action_strings.size());
	print_result_float("Original size", size_original);

	// apply Boolean reductions
	apply_binary_rules(input_br);
	float size_2a = get_size_round(input_br.size(), NSDI_BOOL_SIZE);
	print_result("1-stage Boolean rules", input_br.size());
	print_result_float("1-stage Boolean size", size_2a);

	OneModeResult omr_3a = process_one_mode(input_br, 2, -1);
	print_OMR("1-stage RX filter, TX fpcheck", omr_3a);

	OneModeResult omr_4a = process_one_mode(input_br, 1, -1);
	print_OMR("1-stage RX action, TX none", omr_4a);

	if (ingress_fname != "") {
		vector<NSDIRule> input_br_ingress;
		boost::unordered_map< string, uint > actions_ingress;
		vector<string> action_strings_ingress;
		read_nsdi_file(ingress_fname, input_br_ingress, actions_ingress, action_strings_ingress);
		apply_binary_rules(input_br_ingress);
		float size_ingress_se = get_size_round(input_br_ingress.size(), NSDI_BOOL_SIZE);
		print_result("2-stage ingress orig rules", input_br_ingress.size());
		print_result_float("2-stage ingress orig size", size_ingress_se);

		OneModeResult omr_3b = process_one_mode(input_br_ingress, 1, -1);
		print_OMR("2-stage RX action", omr_3b);
		OneModeResult omr_3c = process_one_mode(input_br_ingress, 0, -1);
		print_OMR("2-stage RX nonconfl", omr_3c);
	}

	if ( (egress0_fname != "") && (egress1_fname != "") ) {
		vector<NSDIRule> input_br_egress_0;
		boost::unordered_map< string, uint > actions_egress_0;
		vector<string> action_strings_egress_0;
		read_nsdi_file(egress0_fname, input_br_egress_0, actions_egress_0, action_strings_egress_0);
		apply_binary_rules(input_br_egress_0);
		OneModeResult omr_egress_4_0 = process_one_mode(input_br_egress_0, 1, -1);
		print_OMR("Divided egress 0", omr_egress_4_0);

		vector<NSDIRule> input_br_egress_1;
		boost::unordered_map< string, uint > actions_egress_1;
		vector<string> action_strings_egress_1;
		read_nsdi_file(egress1_fname, input_br_egress_1, actions_egress_1, action_strings_egress_1);
		apply_binary_rules(input_br_egress_1);
		OneModeResult omr_egress_4_1 = process_one_mode(input_br_egress_1, 1, -1);
		print_OMR("Divided egress 0", omr_egress_4_1);
		
		print_result("Total divided egress rules", input_br_egress_0.size() + input_br_egress_1.size());
		print_result_float("Total divided egress size", get_size_round(input_br_egress_0.size() + input_br_egress_1.size(), NSDI_BOOL_SIZE));
	}

	// uint size_2b_egress = get_size_round(input_br_egress_0.size(), NSDI_BOOL_SIZE) + get_size_round(input_br_egress_1.size(), NSDI_BOOL_SIZE);
	// ofs_t4 << " & " << size_ingress_se << " & " << size_2b_egress << " & " << size_2b_egress+size_ingress_se;
	// ofs_t4 << " & " << pair_3a.first << " & " << pair_3a.second << " & " << pair_3a.first + pair_3a.second;
	// ofs_t4 << " & " << pair_3b.first << " & " << pair_3b.second << " & " << pair_3b.first + pair_3b.second;
	// ofs_t4 << " & " << pair_3c.first << " & " << pair_3c.second << " & " << pair_3c.first + pair_3c.second;
	// ofs_t4 << " & " << pair_4a.first;
	// ofs_t4 << " & " << pair_4b.first << " & " << pair_egress_4_0.first + pair_egress_4_1.first << " & " << pair_egress_4_0.first + pair_egress_4_1.first + pair_4b.first;
	// ofs_t4 << " & " << pair_4c.first << " & " << pair_egress_4_0.first + pair_egress_4_1.first << " & " << pair_egress_4_0.first + pair_egress_4_1.first + pair_4c.first;
}

int main(int argc, char* argv[]) {
	srand(42);

	po::options_description desc("Program options");
	desc.add_options()
	    ("help,?", "produce help message")
	    ("max_width,w", po::value<int>(&max_width)->default_value(-1), "max group width in bits")
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

	process_file_allmodes(in_fname, bits_toremove);
}

