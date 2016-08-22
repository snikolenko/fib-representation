//============================================================================
// Name        : classbench.cpp
// Author      : Sergey Nikolenko
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

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
namespace po = boost::program_options;

const vector<string> fieldnames = { "src IP", "dst IP", "src port", "dst port", "protocol", "???" };

int max_threads;
string input_file, input_file_ranges;

size_t get_share_size(const vector<size_t> v, double share) {
	size_t total = accumulate(v.begin(), v.end(), 0);
	size_t curtotal = 0;
	for (size_t i=0; i<v.size(); ++i) {
		curtotal += v[i];
		if (curtotal >= share * total) return (i+1);
	}
	return v.size();
}

size_t get_thresh_size(const vector<size_t> v, int thresh) {
	size_t i=v.size()-1;
	for (; i>=0; --i) {
		if (v[i] > (uint)thresh) break;
	}
	return (v.size()-i-1);
}

uint field_widths[6] = {32, 32, 16, 16, 8, 16};
uint total_width = 32+32+16+16+8+16;

// uint field_widths[6] = {32, 32, 16, 16, 4, 2};
// uint total_width = 32+32+16+16+4+2;


vector<size_t> get_max_set(const vector<Rule> & v, std::function<bool (const size_t, const size_t)> do_they_intersect) {
	vector<size_t> max_set;
	vector<size_t> have_tried;
	for (size_t i=0; i<v.size(); ) {
		have_tried.push_back(i);
		vector< size_t > set;
		set.push_back(i);
		for (size_t j=0; j<v.size(); ++j) {
			bool ok = all_of(set.begin(), set.end(), [&](size_t h) {
				return !do_they_intersect(h, j);
			});
			if (ok) { 
				set.push_back(j);
			}
		}
		if (set.size() > max_set.size()) {
			max_set = set;
		}
		bool ok = false;
		for (size_t j=i+1; j<v.size(); ++j) {
			ok = all_of(have_tried.begin(), have_tried.end(), [&](size_t h) {
				return do_they_intersect(h, j);
			});
			if (ok) {
				i = j; break;
			}
		}
		// if we got to "default" rules we are out of here
		if (!ok || ( v[i].src_ip_l == 0 && v[i].src_ip_u == 0XFFFF ) || ( v[i].dst_ip_l == 0 && v[i].dst_ip_u == 0xFFFF ) ) break;
	}
	return max_set;
}

bool can_remove_field(const vector<Rule> & v, size_t fld_num, uint max_set_size) {
	vector<size_t> max_set_wo = get_max_set(v, [&](const size_t r1, const size_t r2) {
		for (uint i=0; i<NUM_FIELDS_USED; ++i) {
			if (i==fld_num) continue;
			if (!v[r1].intersects(i, v[r2])) return false;
		}
		return true;
	});
	return (max_set_wo.size() == max_set_size);
}

uint check_and_remove_nonip_fields(vector<Rule> & v, uint max_set_size) {
	uint used_bits = 0;
	for (uint fld_num = 2; fld_num < 6; ++fld_num) {
		if (can_remove_field(v, fld_num, max_set_size)) {
			LOG("Can remove field " << fld_num << " completely.");
			for (uint i=0; i<v.size(); ++i) {
				v[i].set_field_to_zero(fld_num);
			}
		} else {
			LOG("Cannot remove field " << fld_num);
			used_bits += field_widths[fld_num];
		}
	}
	return used_bits;
}

vector<Rule> set_dontcare_in_vector(const vector<Rule> & v, uint iField, uint iBitBegin, uint numBits) {
	vector<Rule> newv;
	for (uint j=0; j<v.size(); ++j) {
		Rule r(v[j], iField, iBitBegin, numBits);
		newv.push_back(r);
 	}
 	return newv;
}

// try to remove as many bits as possible from v starting from iBitBegin for numBits bits
void test_field(vector<Rule> & v, bool **intersects_wrt_not_iField, uint iField, int max_set_size,
		vector<uint> & used_bits_bitwise, size_t iBitBegin, size_t numBits) {
	vector<Rule> newv = set_dontcare_in_vector(v, iField, iBitBegin, numBits);
	vector<size_t> new_max_set = get_max_set(newv, [&](const size_t r1, const size_t r2) {
		return intersects_wrt_not_iField[r1][r2] && newv[r1].intersects(iField, newv[r2]);
	});
	if (new_max_set.size() == max_set_size) {
		LOG("\tcutting bits from " << iBitBegin << " to " << (iBitBegin+numBits-1) << " at once.");
		for (size_t j=0; j<numBits+1; ++j) {
			used_bits_bitwise[j] -= numBits;
		}
		v = newv;
	} else {
		if (numBits == 1) return;
	    test_field(v, intersects_wrt_not_iField, iField, max_set_size, used_bits_bitwise, iBitBegin, (numBits / 2));
	    test_field(v, intersects_wrt_not_iField, iField, max_set_size, used_bits_bitwise, iBitBegin + (numBits / 2), (numBits / 2));
	}
}

void test_field(vector<Rule> & v, bool **intersects_wrt_not_iField, uint iField, int max_set_size, vector<uint> & used_bits_bitwise) {
	test_field(v, intersects_wrt_not_iField, iField, max_set_size, used_bits_bitwise, 0, (iField < 2) ? 32 : 16);
}

vector<size_t> greedyMGR(vector<Rule> & v, const vector<size_t> & subset, int version) {
	for (size_t i=0; i<v.size(); ++i) {
		v[i].spoken_for = true;
	}
	for (size_t j=0; j<subset.size(); ++j) {
		v[subset[j]].spoken_for = false;
	}
	vector<size_t> set_sizes;
	for (size_t i=0; i<v.size(); ++i) {
		if (v[i].spoken_for) continue;
		
		vector< vector< size_t > > sets;
		size_t max_field = 0;
		stringstream ss;

		if (version == 1) {
			sets.resize(NUM_FIELDS_USED);
			
			// one-field greedy algorithm
			for (size_t k=0; k<NUM_FIELDS_USED; ++k) { 				// choose field
				sets[k].push_back(i);
				for (size_t j=i+1; j<v.size(); ++j) {
					if (v[j].spoken_for) continue;
					bool ok = true;
					for (size_t l=0; l<sets[k].size(); ++l) {
						if (v[j].intersects(k, v[sets[k][l]])) {
							ok = false; break;
						}
					}
					if (ok) sets[k].push_back(j);
				}
				if (sets[k].size() > sets[max_field].size()) max_field = k;
			}
			ss << "\twrt field\t" << max_field << " (" << fieldnames[max_field] << ")";
		}

		if (version == 2) {
			sets.resize(NUM_FIELDS_USED * NUM_FIELDS_USED);

			// two-field greedy algorithm
			for (size_t k=0; k<NUM_FIELDS_USED; ++k) {
				for (size_t kk=k+1; kk<NUM_FIELDS_USED; ++kk) { 				// choose two fields
					size_t ind = k * NUM_FIELDS_USED + kk;
					sets[ind].push_back(i);
					for (size_t j=i+1; j<v.size(); ++j) {
						if (v[j].spoken_for) continue;
						bool ok = true;
						for (size_t l=0; l<sets[ind].size(); ++l) {
							if (v[j].intersects(k, kk, v[sets[ind][l]])) {
								ok = false; break;
							}
						}
						if (ok) sets[ind].push_back(j);
					}
					if (sets[ind].size() > sets[max_field].size()) max_field = ind;
				}
			}
			ss << "\twrt fields\t" << (max_field / NUM_FIELDS_USED) << " (" << fieldnames[(max_field / NUM_FIELDS_USED)] << ")\t"
								   << (max_field % NUM_FIELDS_USED) << " (" << fieldnames[(max_field % NUM_FIELDS_USED)] << ")";
		}

		LOG("\tset of size\t" << sets[max_field].size() << "\tstarting from\t" << i << ss.str());
		for (size_t l=0; l<sets[max_field].size(); ++l) {
			v[sets[max_field][l]].spoken_for = true;
		}
		set_sizes.push_back(sets[max_field].size());
	}
	sort(set_sizes.begin(), set_sizes.end(), std::greater<int>());
	return set_sizes;
}

int main(int argc, char* argv[]) {
	srand(time(NULL));

	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,?", "produce help message")
	    ("max_threads,t", po::value<int>(&max_threads)->default_value(1), "number of threads to run")
	    ("input,i", po::value<string>(&input_file)->default_value("/home/snikolenko/networking/ClassBench/test.txt"), "input file")
	    ("addranges,r", po::value<string>(&input_file_ranges)->default_value(""), "input file with additional ranges")
	    ("mark_maxois,m", "mark max OI set and exit")
	    ("bitfields,b", "use bitwise fields and minimize reasonable width")
	    ("cisco,c", "use Cisco input format (not really, but our conversion to pseudo-classbench)")
	    ("mindnf,d", "use Boolean MinDNF techniques on the OIS")
	    ("minac0,a", "use Boolean minimization techniques on both OIS and order-dependend part")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << desc << "\n";
	    return 1;
	}

	bool cisco = vm.count("cisco");

	vector<Rule> v;

	LOG("Reading file " << input_file);
	vector<string> tok;
	vector<size_t> max_set;
	process_file_buffered(input_file, [&] (size_t i, const string & line) {
		boost::algorithm::split(tok, line, boost::is_any_of("\t"), boost::token_compress_on);
		if (tok.size() < NUM_FIELDS) return;
		Rule r(tok, cisco);
		if (tok.size() > NUM_FIELDS) {
			if (tok[NUM_FIELDS][0] == '+') max_set.push_back(i);
		}
		v.push_back(r);
	});
	LOG("Read " << v.size() << " rules. ");

	bool **intersects_wrt_01 = allocate_2d_with_default<bool>(v.size(), v.size(), false);
	for (size_t i=0; i < v.size(); ++i) {
		intersects_wrt_01[i][i] = true;
		for (size_t j=i+1; j < v.size(); ++j) {
			if (v[i].intersects(0, 1, v[j])) {
				intersects_wrt_01[i][j] = true;
				intersects_wrt_01[j][i] = true;
			}
		}
	}
	LOG("Filled indep. w.r.t. fields {0,1} array. Finding max order independent set with greedy algorithm.");
	
	if (!max_set.size()) {
		max_set = get_max_set(v, [&](const size_t r1, const size_t r2) {
			return intersects_wrt_01[r1][r2] && v[r1].intersects(2, 3, v[r2]) && v[r1].intersects(4, 5, v[r2]);
		});
	}

	LOG("Found max order independent set of size " << max_set.size());

	if (vm.count("bitfields")) { // greedily cut out parts of each field
		vector<Rule> v_tmp = v;

		// we begin with non-IP fields, trying to remove each of them completely
		uint used_bits = check_and_remove_nonip_fields(v_tmp, max_set.size());
		uint used_bits_fieldwise = used_bits;

		// if (!can_remove_field(v_tmp, 0, max_set.size())) {
		// 	used_bits_fieldwise += 32;
		// }
		// if (!can_remove_field(v_tmp, 1, max_set.size())) {
		// 	used_bits_fieldwise += 32;
		// }

		vector<uint> used_bits_bitwise = vector<uint>(17, used_bits);
		for (uint iField=0; iField<2; ++iField) {
			LOG("Working with field " << iField);
			used_bits += field_widths[iField];
			for (uint i=0; i<used_bits_bitwise.size(); ++i) { used_bits_bitwise[i] += field_widths[iField]; }

			bool **intersects_wrt_not_iField = allocate_2d_with_default<bool>(v_tmp.size(), v_tmp.size(), true);
			for (size_t i=0; i < v_tmp.size(); ++i) {
				intersects_wrt_not_iField[i][i] = true;
				for (size_t j=i+1; j < v_tmp.size(); ++j) {
					for (size_t k=0; k<NUM_FIELDS_USED; ++k) {
						if ( (k!=iField) && (!v_tmp[i].intersects(k, v_tmp[j])) ) {
							intersects_wrt_not_iField[i][j] = false;
							intersects_wrt_not_iField[j][i] = false;
							break;
						}
					}
				}
			}
			LOG("\tarray ready");

			test_field(v_tmp, intersects_wrt_not_iField, iField, max_set.size(), used_bits_bitwise);

			// bool changed = true;
			// uint iBitSize = (iField < 2) ? 32 : 16;
			// vector<bool> bit_mask(iBitSize, true);
			// while (changed) {
			// 	changed = false;
			// 	int max_new_setsize = -1, max_new_ind = -1;
			// 	for (uint iBit=0; iBit < iBitSize; ++iBit) {
			// 		if (!bit_mask[iBit]) continue;
			// 		vector<Rule> newv;
			// 		for (uint j=0; j<v_tmp.size(); ++j) {
			// 			Rule r(v_tmp[j], iField, iBit);
			// 			newv.push_back(r);
			// 	 	}
			// 		vector<size_t> new_max_set = get_max_set(newv, [&](const size_t r1, const size_t r2) {
			// 			return intersects_wrt_not_iField[r1][r2] && newv[r1].intersects(iField, newv[r2]);
			// 		});
			// 		if ((int)new_max_set.size() > max_new_setsize) {
			// 			max_new_setsize = new_max_set.size();
			// 			max_new_ind = iBit;
			// 			if (max_new_setsize == (int)max_set.size()) break;
			// 		}
			// 	}
			// 	if (max_new_setsize == (int)max_set.size()) {
			// 		changed = true;
			// 		LOG("\tcutting bit " << max_new_ind << " from field " << iField << ", max set of size " << max_new_setsize);
			// 		bit_mask[max_new_ind] = false;
			// 		for (uint j=0; j<v_tmp.size(); ++j) {
			// 			v_tmp[j] = Rule(v_tmp[j], iField, max_new_ind);
			// 		}
			// 		used_bits--;
 		// 			} else {
			// 		LOG("\tbest bit " << max_new_ind << " produces " << max_new_setsize << ", alas.");
			// 	}
			// }

			delete_2d<bool>(v.size(), intersects_wrt_not_iField);
		}

		LOG("Final count: " << used_bits << " bits out of " << total_width);
		LOG("\t" << used_bits_bitwise[0] << " " << used_bits_bitwise[2] << " " << used_bits_bitwise[4] << " " << used_bits_bitwise[8] << " " << used_bits_bitwise[16] << " " << used_bits_fieldwise);
		cout << input_file << " & " << used_bits_bitwise[0] << " " << used_bits_bitwise[2] << " " << used_bits_bitwise[4] << " " << used_bits_bitwise[8] << " " << used_bits_bitwise[16] << " " << used_bits_fieldwise << endl;
		exit(0);
	}

	if (vm.count("mark_maxois")) {
		LOG("Marking MaxOIS in file " << input_file << ".oisnoranges.txt...");
		size_t set_index = 0;
		ofstream ofs(input_file + ".oisnoranges.txt");
		process_file_buffered(input_file, [&] (size_t i, const string & line) {
			if (max_set[set_index] == i) {
				ofs << line << "\t+" << "\n";
				++set_index;
			} else {
				ofs << line << "\t-" << "\n";
			}
		});
		ofs.close();
		LOG("All done.");
		exit(0);
	}

	if (vm.count("mindnf")) {
		bool use_ac0 = vm.count("minac0");
		LOG( (use_ac0 ? "Solving MiNAC0 for the set of rules. Preparing labeled DNF..." : "Solving MinDNF for MaxOIS. Preparing DNF...") );
		vector<BooleanRule *> br;
		size_t set_index = 0;
		process_csv_file_buffered(input_file, "\t", [&] (size_t i, const vector<string> & tok) {
		 	br.push_back(new BooleanRule(tok, cisco));
		});

		if (input_file_ranges != "") {
			LOG("We have additional ranges. Had " << br.size() << " rules before adding ranges...");
			vector<BooleanRule *> br_new;
			set_index = 0;
			vector<string> tok1, tok2;
			process_csv_file_buffered(input_file_ranges, "\t", [&] (size_t iBr, const vector<string> & tok) {
				bool is_in_ois = false;
				if (max_set[set_index] == iBr) {
					is_in_ois = true;
					++set_index;
				} else {
					if (!use_ac0) return;
				}
				// LOG(iBr << ": " << br[iBr]->print());
				boost::split(tok1,tok[2],boost::algorithm::is_any_of(", []'"),boost::token_compress_on);
				boost::split(tok2,tok[3],boost::algorithm::is_any_of(", []'"),boost::token_compress_on);
				for (uint i=0; i<tok1.size(); ++i) {
					if (tok1[i].size() == 0) continue;
					for (uint j=0; j<tok2.size(); ++j) {
						if (tok2[j].size() == 0) continue;
						// LOG(tok1[i] << "\t" << tok2[j]);
						BooleanRule *brule = new BooleanRule(*br[iBr], 88, tok1[i], 104, tok2[j], is_in_ois);
						// LOG(setw(3) << br_new.size() << ": " << brule->print());
						br_new.push_back(brule);
					}
				}
			});
			br = br_new;
		}

		LOG("Starting MinDNF/MinAC0 heuristics for " << br.size() << " rules.");
		uint original_rules_num = br.size();
		bool changed = true;
		int to_be_resolved = -1;
		vector<bool> to_remove(br.size(), false);
		vector<BooleanRule *> to_add;
		size_t iter_count = 0;
		size_t previous_br_size = 0;
		while (changed) {
			LOG("MinDNF iteration " << iter_count << "...");
			++iter_count;
			changed = false;
			for (size_t i=0; i < br.size(); ++i) {
				if (to_remove[i]) continue;
				
				const BooleanRule & r = *br[i];
				for (size_t j=max(previous_br_size, i+1); j < br.size(); ++j) {
					if (to_remove[j]) continue;
					const BooleanRule & s = *br[j];

					// subsumption rule
					switch (r.subsumes(s)) {
						case 1:
							LOG("\trule " << i << " subsumes rule " << j);
							LOG("" << r.print());
							LOG("" << s.print());
							to_remove[j] = true;
							changed = true;
							break;
						case -1:
							LOG("\trule " << j << " subsumes rule " << i);
							LOG("" << r.print());
							LOG("" << s.print());
							to_remove[i] = true;
							changed = true;
							break;
						default: break;
					}

					// resolution rule
					to_be_resolved = r.can_resolve_by(s);
					if (to_be_resolved >= 0) {
						LOG("\trules " << i << " and " << j << " can be resolved:");
						LOG("" << r.print());
						LOG("" << s.print());
						to_add.push_back(new BooleanRule(r, to_be_resolved));
						to_remove[i] = true;
						to_remove[j] = true;
						changed = true;
					}

					if (to_remove[i]) break; // no need to go on if we have decided to remove this rule
				}
			}
			if (!changed) break;

			// adding new rules
			previous_br_size = br.size();
			br.insert(br.end(), to_add.begin(), to_add.end());
			to_add.clear();
			to_remove.resize(br.size(), false);
		}

		vector<BooleanRule*> br_out;
		br_out.reserve(br.size());
		for (size_t i=0; i<br.size(); ++i) {
			if (to_remove[i]) continue;
			br_out.push_back(br[i]);
		}
		br.clear();

		vector<bool> result_mask(BOOL_SIZE, false);
		vector< pair<bool, char> > all_the_same;
		for (size_t i=0; i<br_out.size(); ++i) {
			br_out[i]->mark_in_result_mask(result_mask);
			if (i == 0) {
				br_out[i]->init_all_the_same(all_the_same);
			} else {
				br_out[i]->mark_in_all_the_same(all_the_same);
			}
		}
		uint result_width = accumulate(result_mask.begin(), result_mask.end(), 0);
		uint result_width_withthesame = 0;
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			result_width_withthesame += !(all_the_same[i].first);
		}
		LOG("MinDNF done with " << br_out.size() << " rules left. The result width is " << result_width << " out of " << BOOL_SIZE);

		// for (uint i=0; i<br_out.size(); ++i) {
		// 	LOG(br_out[i]->print());
		// }

		cout << original_rules_num << " & " << BOOL_SIZE << " & " << br_out.size() << " & " << result_width << " & " << result_width_withthesame << " \\\\\n";
		LOG("All done");
		exit(0);
	}

	const vector< vector< vector<size_t> > > fieldsets = {
		{
			{ }
		},
		{
			{2},
			{3},
			{4},
			{5}
		},
		{
			{2, 3},
			{2, 4},
			{2, 5},
			{3, 4},
			{3, 5},
			{4, 5}
		},
		{
			{2, 3, 4},
			{3, 4, 5},
			{2, 4, 5},
			{2, 3, 5}
		}
	};
	vector<size_t> min_fields = {0, 1, 2, 3, 4, 5};
	for (size_t f=0; f < fieldsets.size(); ++f) {
		for (size_t g=0; g<fieldsets[f].size(); ++g) {
			bool is_independent = true;
			for (size_t i=0; i < max_set.size(); ++i) {
				for (size_t j=i+1; j < max_set.size(); ++j) {
					if (intersects_wrt_01[max_set[i]][max_set[j]] && v[max_set[i]].intersects(fieldsets[f][g], v[max_set[j]])) {
						is_independent = false;
						break;
					}
				}
				if (!is_independent) break;
			}
			if (is_independent) {
				min_fields = fieldsets[f][g];
				break;
			}
		}
		if (min_fields.size() < NUM_FIELDS_USED - 2) {
			break;
		} else {
			LOG("\tnot " << (f+2));
		}
	}
	stringstream ss_minfields;
	ss_minfields << "0, 1" << (min_fields.size() > 0 ? ", " : "");
	for (uint i=0; i < min_fields.size(); ++i) ss_minfields << min_fields[i] << (i == min_fields.size() - 1 ? "" : ", ");
	LOG("FSM = { " << ss_minfields.str() << " }");

	// СДЕЛАТЬ ONE-FIELD и TWO-FIELD MGR на ORDER-INDEPENDENT ЧАСТИ

	LOG("Finding max indep set w.r.t. first two fields");
	vector<size_t> max_set_2 = get_max_set(v, [&](const size_t r1, const size_t r2) {
		return intersects_wrt_01[r1][r2];
	});
	LOG("Found max order independent set w.r.t. fields {0, 1} of size " << max_set_2.size());

	delete_2d<bool>(v.size(), intersects_wrt_01);

	vector<size_t> all_inclusive_set;
 	for(size_t i = 0; i < v.size(); ++i ) {
   		all_inclusive_set.push_back(i);
 	}

	LOG("Read " << v.size() << " rules. Starting greedy one-field algorithm.");
	vector<size_t> set_sizes = greedyMGR(v, all_inclusive_set, 1);

	LOG("Greedy one-field algorithm finished with " << set_sizes.size() << " sets, with "
		<< get_share_size(set_sizes, 0.95) << " sets covering 95\% and "
		<< get_share_size(set_sizes, 0.99) << " sets covering 99\% of the rules.");

	LOG("\t" << get_thresh_size(set_sizes, 2) << "\tgroups of size <=2,\t" << get_thresh_size(set_sizes, 5) << "\tgroups of size <=5");


	LOG("Starting greedy two-field algorithm.");
	vector<size_t> set_sizes_2 = greedyMGR(v, all_inclusive_set, 2);

	LOG("Greedy two-field algorithm finished with " << set_sizes_2.size() << " sets, with "
		<< get_share_size(set_sizes_2, 0.95) << " sets covering 95\% and "
		<< get_share_size(set_sizes_2, 0.99) << " sets covering 99\% of the rules.");

	LOG("\t" << get_thresh_size(set_sizes_2, 2) << "\tgroups of size <=2,\t" << get_thresh_size(set_sizes_2, 5) << "\tgroups of size <=5");

	LOG("Starting greedy one-field algorithm on MaxOIS set.");
	vector<size_t> set_sizes_maxois = greedyMGR(v, max_set, 1);

	LOG("Greedy one-field algorithm finished with " << set_sizes_maxois.size() << " sets, with "
		<< get_share_size(set_sizes_maxois, 0.95) << " sets covering 95\% and "
		<< get_share_size(set_sizes_maxois, 0.99) << " sets covering 99\% of the rules.");

	LOG("\t" << get_thresh_size(set_sizes_maxois, 2) << "\tgroups of size <=2,\t" << get_thresh_size(set_sizes_maxois, 5) << "\tgroups of size <=5");

	LOG("Starting greedy two-field algorithm on MaxOIS set.");
	vector<size_t> set_sizes_maxois_2 = greedyMGR(v, max_set, 2);

	LOG("Greedy two-field algorithm finished with " << set_sizes_maxois_2.size() << " sets, with "
		<< get_share_size(set_sizes_maxois_2, 0.95) << " sets covering 95\% and "
		<< get_share_size(set_sizes_maxois_2, 0.99) << " sets covering 99\% of the rules.");

	LOG("\t" << get_thresh_size(set_sizes_maxois_2, 2) << "\tgroups of size <=2,\t" << get_thresh_size(set_sizes_maxois_2, 5) << "\tgroups of size <=5");


	cout << input_file << "\t& " << v.size() << "\t& "
		 << max_set.size() <<  "\t& " << ss_minfields.str() << "\t& " << max_set_2.size() << "\t& "
		 << set_sizes.size() << "\t& " << get_share_size(set_sizes, 0.95) << "\t& " << get_share_size(set_sizes, 0.99) << "\t& "
		 << get_thresh_size(set_sizes, 2) << "\t& " << get_thresh_size(set_sizes, 5) << "\t& "
		 << set_sizes_2.size() << "\t& " << get_share_size(set_sizes_2, 0.95) << "\t& " << get_share_size(set_sizes_2, 0.99) << "\t& "
		 << get_thresh_size(set_sizes_2, 2) << "\t& " << get_thresh_size(set_sizes_2, 5) << "\t& "
		 << set_sizes_maxois.size() << "\t& " << get_share_size(set_sizes_maxois, 0.95) << "\t& " << get_share_size(set_sizes_maxois, 0.99) << "\t& "
		 << get_thresh_size(set_sizes_maxois, 2) << "\t& " << get_thresh_size(set_sizes_maxois, 5) << "\t& "
		 << set_sizes_maxois_2.size() << "\t& " << get_share_size(set_sizes_maxois_2, 0.95) << "\t& " << get_share_size(set_sizes_maxois_2, 0.99) << "\t& "
		 << get_thresh_size(set_sizes_maxois_2, 2) << "\t& " << get_thresh_size(set_sizes_maxois_2, 5) << "\\\\" << endl;
}

