#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>


#define NUM_FIELDS 6
#define NUM_FIELDS_USED 6

#define NSDI_NUM_FIELDS 1

#define NSDI_BOOL_SIZE 32
// #define NSDI_BOOL_SIZE 128
// #define RANDOM_PADDING

#define NSDI_BOOL_SIZE_IPV6 128

#define BOOL_SIZE 64
// #define BOOL_SIZE 98

using namespace std;

bool intervals_intersect(uint l1, uint u1, uint l2, uint u2) {
	return (u1 >= l2) && (u2 >= l1);
}

// tmp variables
unsigned char ip[4], cisco_mask[4], cisco_ip_l[4], cisco_ip_u[4], mask, prot_m, prot_l, ipv6_field;
unsigned int field_m, field_v;
vector<string> ipv6_tok;
char ipv6[40];

template<typename T>
T *allocate_1d_with_default(size_t dim, T def) {
    T *r = new T[dim];
    for (size_t i=0; i < dim; ++i) {
        r[i] = def;
    }
    return r;
}

template<typename T>
T **allocate_2d_with_default(size_t dim1, size_t dim2, T def) {
    T **r = new T*[dim1];
    for (size_t i=0; i < dim1; ++i) {
        r[i] = allocate_1d_with_default(dim2, def);
    }
    return r;
}

template<typename T>
void delete_2d(size_t dim_one, T **r) {
    for (size_t i=0; i<dim_one; ++i) {
        delete r[i];
    }
    delete r;
}

void scan_ip_with_mask(const char * src_string, uint * ip_l, uint * ip_u) {
	sscanf(src_string, "%hhu.%hhu.%hhu.%hhu/%hhu", &ip[0], &ip[1], &ip[2], &ip[3], &mask);
	*ip_l = (( (ip[0] << 24) + (ip[1] << 16) + (ip[2] << 8) + ip[3] ) >> (32-mask)) << (32-mask);
	*ip_u = *ip_l + (1 << (32-mask)) - 1;
}

void scan_cisco_ip_with_mask(const char * src_string, uint * ip_l, uint * ip_u) {
	sscanf(src_string, "%hhu.%hhu.%hhu.%hhu/%hhu.%hhu.%hhu.%hhu",
		&ip[0], &ip[1], &ip[2], &ip[3], &cisco_mask[0], &cisco_mask[1], &cisco_mask[2], &cisco_mask[3]);
	for (uint i=0; i<4; ++i) {
		cisco_ip_l[i] = ip[i] & cisco_mask[3-i];
		cisco_ip_u[i] = cisco_ip_l[i] + (!cisco_mask[3-i]);
		// LOG((int)ip[i] << "/" << (int)cisco_mask[3-i] << "\t=\t" << (int)cisco_ip_l[i] << " : " << (int)cisco_ip_u[i]);
	}
	*ip_l = (cisco_ip_l[0] << 24) + (cisco_ip_l[1] << 16) + (cisco_ip_l[2] << 8) + cisco_ip_l[3];
	*ip_u = (cisco_ip_u[0] << 24) + (cisco_ip_u[1] << 16) + (cisco_ip_u[2] << 8) + cisco_ip_u[3];
}

void scan_ip_without_mask_boolean(const char * src_string, char * bitstring) {
	if (src_string[0] == '@') {
		sscanf(src_string, "@%hhu.%hhu.%hhu.%hhu", &ip[0], &ip[1], &ip[2], &ip[3]);
	} else {
		sscanf(src_string, "%hhu.%hhu.%hhu.%hhu", &ip[0], &ip[1], &ip[2], &ip[3]);
	}
	for (size_t i=0; i<4; ++i) {
		for (size_t j=0; j<8; ++j) {
			bitstring[8*i+j] =  ( (ip[i] & (1 << j)) == 0 ? 0 : 1);
		}
	}
}

void scan_ip_with_mask_boolean(const char * src_string, char * bitstring) {
	if (src_string[0] == '@') {
		sscanf(src_string, "@%hhu.%hhu.%hhu.%hhu/%hhu", &ip[0], &ip[1], &ip[2], &ip[3], &mask);
	} else {
		sscanf(src_string, "%hhu.%hhu.%hhu.%hhu/%hhu", &ip[0], &ip[1], &ip[2], &ip[3], &mask);
	}
	for (size_t i=0; i<4; ++i) {
		for (size_t j=0; j<8; ++j) {
			bitstring[8*i+7-j] =  ( (ip[i] & (1 << j)) == 0 ? 0 : 1);
		}
	}
	for (size_t j=mask; j<32; ++j) {
		bitstring[j] = '*';
	}
}

void scan_16bit_number_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%x", &field_v);
	for (size_t j=0; j<16; ++j) {
		bitstring[j] =  ( (prot_l & (1 << j)) == 0 ? 0 : 1);
	}
}

void scan_16bit_ipv6(const char * src_string, char * bitstring) {
	sscanf(src_string, "%x", &field_v);
	for (size_t j=0; j<16; ++j) {
		bitstring[16-j] =  ( (field_v & (1 << j)) == 0 ? 0 : 1);
	}
}

void scan_ipv6_with_mask_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%[0123456789abcdef:]/%hhu", ipv6, &mask);
	boost::algorithm::split(ipv6_tok, src_string, boost::is_any_of(":"));
	uint cur_group = 0;
	for (size_t i=0; i<ipv6_tok.size(); ++i) {
		if (ipv6_tok[i].size() == 0) {
			cur_group = 8 - ipv6_tok.size();
		} else {
			scan_16bit_ipv6(ipv6_tok[i].c_str(), bitstring + 8*cur_group);
			++cur_group;
		}
	}
	for (size_t j=mask; j<128; ++j) {
		bitstring[j] = '*';
	}
}

void scan_cisco_ip_with_mask_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%hhu.%hhu.%hhu.%hhu/%hhu.%hhu.%hhu.%hhu",
		&ip[0], &ip[1], &ip[2], &ip[3], &cisco_mask[0], &cisco_mask[1], &cisco_mask[2], &cisco_mask[3]);
	// LOG((int)ip[0] << "\t" << (int)cisco_mask[3]);
	for (size_t i=0; i<4; ++i) {
		for (size_t j=0; j<8; ++j) {
			bitstring[8*i+j] =  ( (ip[i] & (1 << j)) == 0 ? 0 : 1);
		}
		for (size_t j=0; j<8; ++j) {
			if ( !(cisco_mask[3-i] & (1 << j)) ) bitstring[8*i+j] = '*';
		}
	}
}

string binary_string(uint x, uint sz=32) {
	ostringstream ss;
	for (int i = sz-1; i >= 0; i--)
    	ss << ((x >> i) & 1);
    return ss.str();
}

uint remove_bit_uint(uint x, size_t b) {
	// LOG("x=" << x << " minus " << b << " = " << ( x - ( ((x >> b) & 1) << b ) ) << "\t" << binary_string(x) << "\t" << binary_string(( x - ( ((x >> b) & 1) << b ) )));
	return ( x - ( ((x >> b) & 1) << b ) );
}

void scan_8bit_field_with_mask_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%hhx/%hhx", &prot_l, &prot_m);
	for (size_t j=0; j<8; ++j) {
		bitstring[j] =  ( (prot_l & (1 << j)) == 0 ? 0 : 1);
	}
	for (size_t j=0; j<8; ++j) {
		if ((prot_m & (1 << j)) == 0) bitstring[j] = '*';
	}
}

void scan_8bit_number_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%hhx", &prot_l);
	for (size_t j=0; j<8; ++j) {
		bitstring[j] =  ( (prot_l & (1 << j)) == 0 ? 0 : 1);
	}
}

void scan_16bit_field_with_mask_boolean(const char * src_string, char * bitstring) {
	sscanf(src_string, "%x/%x", &field_v, &field_m);
	for (size_t j=0; j<16; ++j) {
		bitstring[j] =  ( (field_v & (1 << j)) == 0 ? 0 : 1);
	}
	for (size_t j=0; j<16; ++j) {
		if ((field_m & (1 << j)) == 0) bitstring[j] = '*';
	}
}

class Rule {
public:
	uint src_ip_l, src_ip_u, dst_ip_l, dst_ip_u;
	uint src_port_l, src_port_u, dst_port_l, dst_port_u;
	unsigned char prot_l, prot_u;
	unsigned int field, field_m;
	bool cisco;
	bool nsdi;
	bool spoken_for;

	Rule() : cisco(false), spoken_for(false) {	}

	Rule(const vector<string> & tok, int inp_format) : cisco(inp_format == 1), nsdi(inp_format == 2), spoken_for(false) {
		if (cisco) {
			scan_cisco_ip_with_mask(tok[0].c_str(), &src_ip_l, &src_ip_u);
			scan_cisco_ip_with_mask(tok[1].c_str(), &dst_ip_l, &dst_ip_u);
		} else {
			scan_ip_with_mask(tok[0].c_str() + 1, &src_ip_l, &src_ip_u);
			scan_ip_with_mask(tok[1].c_str(), &dst_ip_l, &dst_ip_u);
		}
		sscanf(tok[2].c_str(), "%d : %d", &src_port_l, &src_port_u);
		sscanf(tok[3].c_str(), "%d : %d", &dst_port_l, &dst_port_u);
		if (cisco) {
			sscanf(tok[4].c_str(), "%hhx", &prot_l);
			sscanf(tok[5].c_str(), "%hhx", &prot_u);
		} else {
			sscanf(tok[4].c_str(), "%hhx/%hhx", &prot_l, &prot_m);
			if (prot_m == (unsigned char)0xFF) {
				prot_u = prot_l;
			} else {					// this means it is 0x00
				prot_l = 0x00;
				prot_u = 0xFF;
			}
			sscanf(tok[5].c_str(), "%u/%u", &field, &field_m);
		}
	}

	Rule(const Rule & r, size_t fInd, size_t bInd) : src_ip_l(r.src_ip_l), src_ip_u(r.src_ip_u), dst_ip_l(r.dst_ip_l), dst_ip_u(r.dst_ip_u),
			src_port_l(r.src_port_l), src_port_u(r.src_port_u), dst_port_l(r.dst_port_l), dst_port_u(r.dst_port_u),
			prot_l(r.prot_l), prot_u(r.prot_u), field(r.field), field_m(r.field_m), cisco(r.cisco), spoken_for(false) {
		// remove one bit from one field
		switch (fInd) {
			case 0: src_ip_l = remove_bit_uint(src_ip_l, bInd); src_ip_u = remove_bit_uint(src_ip_u, bInd); break;
			case 1: dst_ip_l = remove_bit_uint(dst_ip_l, bInd); dst_ip_u = remove_bit_uint(dst_ip_u, bInd); break;
			case 2: src_port_l = remove_bit_uint(src_port_l, bInd); src_port_u = remove_bit_uint(src_port_u, bInd); break;
			case 3: dst_port_l = remove_bit_uint(dst_port_l, bInd); dst_port_u = remove_bit_uint(dst_port_u, bInd); break;
		}
	}

	Rule(const Rule & r, size_t fInd, size_t bInd, size_t sz) :
			src_ip_l(r.src_ip_l), src_ip_u(r.src_ip_u), dst_ip_l(r.dst_ip_l), dst_ip_u(r.dst_ip_u),
			src_port_l(r.src_port_l), src_port_u(r.src_port_u), dst_port_l(r.dst_port_l), dst_port_u(r.dst_port_u),
			prot_l(r.prot_l), prot_u(r.prot_u), field(r.field), field_m(r.field_m), cisco(r.cisco), spoken_for(false) {
		// remove several bits from one field
		for (size_t i = bInd; i < bInd + sz; ++i) {
			switch (fInd) {
				case 0:	src_ip_l = remove_bit_uint(src_ip_l, i); src_ip_u = remove_bit_uint(src_ip_u, i); break;
				case 1: dst_ip_l = remove_bit_uint(dst_ip_l, i); dst_ip_u = remove_bit_uint(dst_ip_u, i); break;
				case 2: src_port_l = remove_bit_uint(src_port_l, i); src_port_u = remove_bit_uint(src_port_u, i); break;
				case 3: dst_port_l = remove_bit_uint(dst_port_l, i); dst_port_u = remove_bit_uint(dst_port_u, i); break;
			}
		}
	}

	bool intersects(size_t k, const Rule & r) const {
		uint m;
		switch (k) {
			case 0: return intervals_intersect(src_ip_l, src_ip_u, r.src_ip_l, r.src_ip_u);
			case 1: return intervals_intersect(dst_ip_l, dst_ip_u, r.dst_ip_l, r.dst_ip_u);
			case 2: return intervals_intersect(src_port_l, src_port_u, r.src_port_l, r.src_port_u);
			case 3: return intervals_intersect(dst_port_l, dst_port_u, r.dst_port_l, r.dst_port_u);
			case 4: return ( cisco ? (prot_l == r.prot_l) : intervals_intersect(prot_l, prot_u, r.prot_l, r.prot_u) );
			case 5:
				if (cisco) {
					return (prot_u == r.prot_u);
				} else {
					m = field_m & r.field_m;
					if (m == 0) return true;
					return (field & m) == (r.field & m);
				}
			default: return false;
		}
	}

	bool intersects(size_t k, size_t kk, const Rule & r) const {
		return intersects(k, r) && intersects(kk, r);
	}

	bool intersects(const Rule & r) const {
		return intersects(0, r) && intersects(1, r) && intersects(2, r) && intersects(3, r) && intersects(4, r) && intersects(5, r);
	}

	bool intersects(const vector<size_t> & f, const Rule & r) const {
		for (size_t i=0; i<f.size(); ++i) {
			if (!intersects(f[i], r)) return false;
		}
		return true;
	}

	void set_field_to_zero(size_t k) {
		switch (k) {
			case 0: src_ip_l = 0; src_ip_u = 0; break;
			case 1: dst_ip_l = 0; dst_ip_u = 0; break;
			case 2: src_port_l = 0; src_port_u = 0; break;
			case 3: dst_port_l = 0; dst_port_u = 0; break;
			case 4: prot_l = 0; prot_u = 0; break;
			case 5: field = 0; field_m = 0; break;
		}
	}

	string print_field(size_t k) const {
		ostringstream ss;
		switch (k) {
			case 0: ss << "[ " << src_ip_l << " : " << src_ip_u << " ]"; break;
			case 1: ss << "[ " << dst_ip_l << " : " << dst_ip_u << " ]"; break;
			case 2: ss << "[ " << src_port_l << " : " << src_port_u << " ]"; break;
			case 3: ss << "[ " << dst_port_l << " : " << dst_port_u << " ]"; break;
			case 4: ss << "[ " << prot_l << " : " << prot_u << " ]"; break;
			case 5: ss << "[ " << field << " / " << field_m << " ]"; break;
		}
		return ss.str();
	}
};


class BooleanRule {
protected:
	char *x;
	bool unit;
	bool is_in_ois;

	void update_vars() {
		unit = false;
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			if (x[i] != '*') {
				if (unit) {
					unit = false; break;
				} else {
					unit = true;
				}
			}
		}
	}
public:
	BooleanRule(const vector<string> & tok, int in_format) : unit(false), is_in_ois(true) {
		x = new char[BOOL_SIZE];
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			x[i] = '*';
		}
		if (in_format == 1) {
			scan_cisco_ip_with_mask_boolean(tok[0].c_str(), x);
			scan_cisco_ip_with_mask_boolean(tok[1].c_str(), x+32);
			scan_8bit_number_boolean(tok[4].c_str(), x+64);
			scan_16bit_number_boolean(tok[5].c_str(), x+72);
		} else if (in_format == 2) {
			scan_ip_with_mask_boolean(tok[0].c_str(), x);
			scan_ip_without_mask_boolean(tok[1].c_str(), x+32);
		} else {
			// LOG(tok[4] << "\t" << tok[5]);
			scan_ip_with_mask_boolean(tok[0].c_str(), x);
			scan_ip_with_mask_boolean(tok[1].c_str(), x+32);
			scan_8bit_field_with_mask_boolean(tok[4].c_str(), x+64);
			scan_16bit_field_with_mask_boolean(tok[5].c_str(), x+72);
			// LOG(print());
		}
		update_vars();
	}

	// construct BooleanRule as a result of resolution
	BooleanRule(const BooleanRule & b, int to_resolve) {
		x = new char[BOOL_SIZE];
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			x[i] = b.x[i];
		}
		if (to_resolve >= 0 && to_resolve < BOOL_SIZE) {
			x[to_resolve] = '*';
		}
	}

	// construct BooleanRule with additional ranges and ois mark
	BooleanRule(const BooleanRule & b, size_t start1, const string & range1, size_t start2, const string & range2, bool ois) : is_in_ois(ois) {
		x = new char[BOOL_SIZE];
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			x[i] = b.x[i];
		}
		for (size_t i=0; i<range1.size(); ++i) {
			x[start1+i] = (range1[i] == '0') ? 0 : ( range1[i]  == '1' ? 1 : '*');
		}
		for (size_t i=0; i<range2.size(); ++i) {
			x[start2+i] = (range2[i] == '0') ? 0 : ( range2[i]  == '1' ? 1 : '*');
		}
	}

	string print() const {
		string res;
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			res += ( x[i] == 0 ? '0' : ( x[i] == 1 ? '1' : '*') );
		}
		return res;
	}

	~BooleanRule() {
		delete [] x;
	}

	// return 0 if none subsumes the other, 1 if this subsumes b, -1 if b subsumes this
	int subsumes(const BooleanRule & b) const {
		int res = -1;
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			if (x[i] != b.x[i]) {
				if (x[i] == '*') {
					if (res >= 0) res = 1;
					else return false;
				}
				if (b.x[i] == '*') {
					if (res <= 0) res = -1;
					else return false;
				}
				return false;
			}
		}
		return res;
	}

	int can_resolve_by(const BooleanRule & b) const {
		int resolve_by = -1;
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			if (x[i] != b.x[i]) {
				if (resolve_by == -1 && x[i] != '*' && b.x[i] != '*') resolve_by = i;
				else return -1;
			}
		}
		return resolve_by;
	}

	void mark_in_result_mask(vector<bool> & result_mask) {
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			if (x[i] != '*') result_mask[i] = true;
		}
	}

	void init_all_the_same(vector< pair<bool, char> > & all_the_same) {
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			all_the_same.push_back(make_pair(true, x[i]));
		}
	}

	void mark_in_all_the_same(vector< pair<bool, char> > & all_the_same) {
		for (size_t i=0; i<BOOL_SIZE; ++i) {
			if (!all_the_same[i].first) continue;
			if (all_the_same[i].second != x[i]) all_the_same[i].first = false;
		}
	}
};


class NSDIRule {
public:
	char x[NSDI_BOOL_SIZE];
	bool unit;
	bool is_in_ois;
	uint a;
	uint l;

	NSDIRule(const string & s, uint act, bool read_binary=false) : unit(false), is_in_ois(true), a(act) {
		if (read_binary) {
			for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
				x[i] = (s[i] == '0' ? 0 : (s[i] == '1' ? 1 : '*') );
			}
		} else {
			if (NSDI_BOOL_SIZE == 128) {
				for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
					x[i] = 0;
				}
	#ifdef RANDOM_PADDING
				for (size_t i=0; i<NSDI_BOOL_SIZE-32; ++i) {
					x[i] = rand() % 2;
				}
				scan_ip_with_mask_boolean(s.c_str(), x + NSDI_BOOL_SIZE - 32);
	#else
				scan_ipv6_with_mask_boolean(s.c_str(), x);
	#endif
			} else {
				for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
					x[i] = '*';
				}
				if (s.size() == NSDI_BOOL_SIZE) {
					for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
						x[i] = s[i] == '0' ? 0 : ( s[i] == '1' ? 1 : '*');
					}
				} else {
					scan_ip_with_mask_boolean(s.c_str(), x);
				}
			}
		}
	}

	NSDIRule(const vector<string> & tok, int in_format) : unit(false), is_in_ois(true) {
		if (NSDI_BOOL_SIZE == 128) {
			for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
				x[i] = 0;
			}
			scan_ipv6_with_mask_boolean(tok[0].c_str(), x);
			cout << tok[0] << this->print() << endl;
		} else {
			for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
				x[i] = '*';
			}
			scan_ip_with_mask_boolean(tok[0].c_str(), x);
		}
	}

	// construct NSDIRule as a result of resolution
	NSDIRule(const NSDIRule & b, int to_resolve) {
		// x = new char[NSDI_BOOL_SIZE];
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			x[i] = b.x[i];
		}
		if (to_resolve >= 0 && to_resolve < NSDI_BOOL_SIZE) {
			x[to_resolve] = '*';
		}
		a = b.a;
		l = (to_resolve == (int)b.l-1) ? b.l - 1 : b.l;
	}

	string print() const {
		string res;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			res += ( x[i] == 0 ? '0' : ( x[i] == 1 ? '1' : '*') );
		}
		return res;
	}

	string print(const vector<bool> & mask) const {
		string res;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			res +=  mask[i] ? ( x[i] == 0 ? '0' : ( x[i] == 1 ? '1' : '*') ) : '.';
		}
		return res;
	}

	~NSDIRule() {
		// delete [] x;
	}

	// return 0 if none subsumes the other, 1 if this subsumes b, -1 if b subsumes this
	int subsumes(const NSDIRule & b) const {
		int res = -1;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if (x[i] != b.x[i]) {
				if (x[i] == '*') {
					if (res >= 0) res = 1;
					else return false;
				}
				if (b.x[i] == '*') {
					if (res <= 0) res = -1;
					else return false;
				}
				return false;
			}
		}
		return res;
	}

	int can_resolve_by(const NSDIRule & b) const {
		if (a != b.a) return -1;
		int resolve_by = -1;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if (x[i] != b.x[i]) {
				if (resolve_by == -1 && x[i] != '*' && b.x[i] != '*') resolve_by = i;
				else return -1;
			}
		}
		return resolve_by;
	}

	void mark_in_result_mask(vector<bool> & result_mask) {
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if (x[i] != '*') result_mask[i] = true;
		}
	}

	void init_all_the_same(vector< pair<bool, char> > & all_the_same) {
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			all_the_same.push_back(make_pair(true, x[i]));
		}
	}

	void mark_in_all_the_same(vector< pair<bool, char> > & all_the_same) {
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if (!all_the_same[i].first) continue;
			if (all_the_same[i].second != x[i]) all_the_same[i].first = false;
		}
	}

	bool intersects(const NSDIRule & b) const {
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if ((x[i] == 0 && b.x[i] == 1) || (x[i] == 1 && b.x[i] == 0)) return false;
		}
		return true;
	}

	// modes:
	// 0 -- relaxed SE
	// 1 -- action order independence
	// 2 -- filter order independence
	void mark_relaxed_intersections(const NSDIRule & b, vector<bool> & mask, uint mode=0) const {
		if ( (mode == 0 || mode == 1) && a == b.a) return;
		int first_difference = -1;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			// if (!mask[i]) continue;
			if ( ((mode == 0) && (( (b.x[i] == 1) && ( x[i] == 0 || x[i] == '*') ) || ( (b.x[i] == 0) && ( x[i] == 1 || x[i] == '*') )))
			|| ( (mode > 0) && ((x[i] == 0 && b.x[i] == 1) || (x[i] == 1 && b.x[i] == 0)) )) {
				if (first_difference >= 0) return;
				else first_difference = i;
			}
		}
		if (first_difference >= 0 && mask[first_difference]) {
			LOG("bit " << first_difference << ":\t" << print() << "\t" << b.print());
			mask[first_difference] = false;
		}
	}

	// modes:
	// 0 -- relaxed SE
	// 1 -- action order independence
	// 2 -- filter order independence
	bool intersects_in_mode_masked_above(const NSDIRule & b, const vector<bool> & mask, uint mode=0) const {
		if ( (mode == 0 || mode == 1) && a == b.a) return false;
		for (size_t i=0; i<NSDI_BOOL_SIZE; ++i) {
			if (!mask[i]) continue;
			if ( ((mode == 0) && (( (b.x[i] == 1) && ( x[i] == 0 || x[i] == '*') ) || ( (b.x[i] == 0) && ( x[i] == 1 || x[i] == '*') )))
			|| ( (mode > 0) && ((x[i] == 0 && b.x[i] == 1) || (x[i] == 1 && b.x[i] == 0)) )) {
				return false;
			}
		}
		return true;
	}

	void mark_relaxed_intersections_masked(const NSDIRule & b,
			vector<bool> & mask,
			const vector<bool> & removed_bits,
			uint min_active_bit, uint max_active_bit,
			uint mode=0) const {
		if ( (mode == 0 || mode == 1) && a == b.a) return;
		int first_difference = -1;
		for (size_t i=min_active_bit; i<max_active_bit; ++i) {
			if (removed_bits[i]) continue;
			if ( ((mode == 0) && (( (b.x[i] == 1) && ( x[i] != 1 ) ) || ( (b.x[i] == 0) && ( x[i] != 0 ) )))
			|| ( (mode > 0) && ((x[i] == 0 && b.x[i] == 1) || (x[i] == 1 && b.x[i] == 0)) )) {
				if (first_difference >= 0) return;
				else first_difference = i;
			}
		}
		if (first_difference >= 0) {
			// LOG("bit " << first_difference << ":\t" << print() << "\t" << b.print());
			mask[first_difference] = true;
		}
	}


};


// import numpy as np
// import glob
// from operator import itemgetter

// for fname in glob.glob('fibs_needed/*.txt'):
// 	print fname
// 	arr = []
// 	mA = {}
// 	with open(fname) as f:
// 	    for line in f:
// 	    	act = line.strip().split()[-1]
// 	    	if not act in mA:
// 	    		mA[act] = len(mA)
// 	        arr.append(mA[act])
// 	arr = np.array(arr)
// 	asum = sorted( [(i, np.sum(arr == i)) for i in xrange(len(mA))], key=itemgetter(1), reverse=True )
// 	asums = [x[1] for x in asum]
// 	msumindex = { asum[i][0] : i for i in xrange(len(asum)) }
// 	for i in xrange(len(asum)):
// 		if np.sum(asums[:i]) <= len(arr) / 2 and np.sum(asums[:i+1]) >= len(arr) / 2:
// 			divider = i+1
// 			break

// 	new_arr = [ 0 if msumindex[i] < divider else 1 for i in arr ]
// 	with open(fname) as f:
// 		with open('fibs_divided/%s' % fname.split('/')[-1], 'w') as outf:
// 			i = 0
// 			for line in f:
// 				cur_rule = line.strip().split()
// 				outf.write('%s %d\n' % (cur_rule[0], new_arr[i]) )
// 				i += 1

// 	new_arr = [ 0 if msumindex[i] < divider else 1 for i in arr ]
// 	with open(fname) as f:
// 		with open('fibs_egress/0.%s' % fname.split('/')[-1], 'w') as outf0:
// 			with open('fibs_egress/1.%s' % fname.split('/')[-1], 'w') as outf1:
// 				i = 0
// 				for line in f:
// 					if new_arr[i] == 0:
// 						outf0.write(line)
// 					else:
// 						outf1.write(line)
// 					i += 1

