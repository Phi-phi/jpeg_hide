#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <typeinfo>
#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <boost/regex.hpp>

#define BYTE unsigned char
#define VECTARR(TYPE) vector< vector<TYPE> >
#define BYTE_ITR vector<BYTE>::iterator
using namespace std;
using namespace boost;

typedef struct dqt {
	int pq;
	int tq;
	vector<BYTE> table;
} dqt_t;

typedef struct haff {
	string haff;
	int r_len;
	int val;
} haff_t;

typedef struct dht {
	int tc;
	int th;
	vector<haff_t> haffs;
} dht_t;

typedef struct portion {
	int id;
	int hn;
	int vn;
	int dqn;
} por_t;

typedef struct sof {
	int p_count;
	int wx;
	int hy;
	vector<por_t> pors;
} sof_t;

typedef struct dh_table_inf {
	int id;
	int act_cd;
	int dct_inf;
} dht_inf;

typedef struct sos {
	int dht_inf_count;
	int s_start;
	int s_end;
	int ah;
	int al;
	vector<dht_inf> infs;
} sos_t;

void indicate(vector<BYTE> vect) {
	for (int i = 0; i < vect.size(); ++i) {
		printf("%2x", vect[i]);
	}
	printf("\n");
}

bool cmp_bin(BYTE vect, string bin) {
	//printf("%d %d\n", vect, strtol(bin.c_str(), NULL, 16));
	return (int)vect == strtol(bin.c_str(), NULL, 16);
}

int search_vct(vector<BYTE> vect, string query) {
	if (query.size() > 2) {
		string first_q = query.substr(0, 2);
		string second_q = query.substr(2);

		for (int i = 0; i < vect.size() - 1; ++i) {
			if (cmp_bin(vect[i], first_q) && cmp_bin(vect[i + 1], second_q)) {
				return i;
			}
		}
		return -1;
	}
	for (int i = 0; i < vect.size(); ++i) {
		if ((char*)&vect[i] == query) {
			return i;
		}
	}
	return -1;
}

VECTARR(BYTE) get_segments(vector<BYTE> vect, string name) {
	int index;
	VECTARR(BYTE) arr;
	while ((index = search_vct(vect, name)) > -1) {
		vector<BYTE> segment;
		int len = (vect[index + 2] << 8) + vect[index + 3] + 2;
		for (int i = index + 2; i < index + len; ++i) {
			segment.push_back(vect[i]);
		}
		arr.push_back(segment);
		vect.erase(vect.begin() + index, vect.begin() + index + len);
	}
	return arr;
}


vector<dqt_t> ana_dqt(VECTARR(BYTE) vect) {
	vector<dqt_t> tables;
	for (int i = 0; i < vect.size(); ++i) {
		dqt_t t;
		vector<BYTE> seg = vect[i];
		int len = (seg[0] << 8) + seg[1] - 2;
		seg.erase(seg.begin(), seg.begin() + 2);
		for (int j = 0; j < len / 65; ++j) {
			t.pq = seg[0] >> 4;
			t.tq = seg[0] & 0x0f;
			for (int k = 1; k < 65; ++k) {
				t.table.push_back(seg[k]);
			}
			seg.erase(seg.begin(), seg.begin() + 65);
		}
		tables.push_back(t);
	}
	return tables;
}

vector<dht_t> ana_dht(VECTARR(BYTE) vect) {
	vector<dht_t> tables;
	for (int i = 0; i < vect.size(); ++i) {
		int index = 0;
		dht_t t;
		vector<BYTE> seg = vect[i];
		int len = (seg[0] << 8) + seg[1] - 2;
		printf("dht len->%d-now %d\n", len, index);
		if (len < 16) {
			printf("Length error.\n shorter than 16.\n");
			exit(-1);
		}
		indicate(seg);
		seg.erase(seg.begin(), seg.begin() + 2);
		indicate(seg);
		int now_haff = 0;
		vector<haff_t> haff_arr;
		t.tc = seg[index] >> 4;
		t.th = seg[index] & 0x0f;
		++index;
		if (!t.tc) printf("DC table Cd:%d\n", t.th);
		else printf("AC table Cd:%d\n", t.th);

		for (int bit_length = 1; bit_length <= 16; ++bit_length) {
			int haff_counts = seg[bit_length];
			if (haff_counts > 0) {
				for (int j = 1; j < haff_counts + 1; ++j) {
					haff_t ht;
					unsigned char bin_haff = now_haff;
					dynamic_bitset<> bit_value(bit_length, now_haff);
					to_string(bit_value, ht.haff);

					++now_haff;
					haff_arr.push_back(ht);
				}
			}
			now_haff = now_haff << 1;
			++index;
			if (index > len) break;
		}
		index = 17;
		for (int h_arr_i = 0; h_arr_i < haff_arr.size(); ++h_arr_i) {
			int val = seg[index];
			if (!t.tc) {
				haff_arr[h_arr_i].r_len = 0;
				haff_arr[h_arr_i].val = val;
			}
			else {
				haff_arr[h_arr_i].r_len = val >> 4;
				haff_arr[h_arr_i].val = val & 0x0f;
			}
			++index;
		}
		t.haffs = haff_arr;
		tables.push_back(t);
	}
	return tables;
}

vector<por_t> ana_sof_block(vector<BYTE> vect, int count) {
	vector<por_t> tables;
	int index = 7;
	for (int i = 0; i < count; ++i) {
		por_t t;
		t.id = vect[++index];
		int hv = vect[++index];
		t.hn = hv >> 4;
		t.vn = hv & 0x0f;
		t.dqn = vect[++index];
		tables.push_back(t);
	}
	return tables;
}

vector<sof_t> ana_sof(VECTARR(BYTE) vect) {
	vector<sof_t> tables;
	for (int i = 0; i < vect.size(); ++i) {
		vector<BYTE> seg = vect[i];
		sof_t t;
		int len = (seg[0] << 8) + seg[1] - 2;
		t.hy = (seg[3] << 8) + seg[4];
		t.wx = (seg[5] << 8) + seg[6];
		t.p_count = seg[7];
		t.pors = ana_sof_block(seg, t.p_count);
		tables.push_back(t);
	}
	return tables;
}

vector<sos_t> ana_sos(VECTARR(BYTE) vect) {
}

int main(int argc, char* argv[]) {
	const char *filepath = argv[1];
	vector<BYTE> img_raw;
	int zig_zag[] = { 0, 1, 8, 16, 9, 2, 3, 10,
		17, 24, 32, 25, 18, 11, 4, 5,
		12, 19, 26, 33, 40, 48, 41, 34,
		27, 20, 13, 6, 7, 14, 21, 28,
		35, 42, 49, 56, 57, 50, 43, 36,
		29, 22, 15, 23, 30, 37, 44, 51,
		58, 59, 52, 45, 38, 31, 39, 46,
		53, 60, 61, 54, 47, 55, 62, 63 };

	ifstream ifs(filepath, ios::in | ios::binary);

	if (!ifs) {
		cout << "FIle open error.\n";
		return 1;
	}

	img_raw.resize(ifs.seekg(0, ios::end).tellg());

	ifs.seekg(0, ios::beg).read((char*)&img_raw[0], img_raw.size());
	ifs.close();


	printf("%d\n", img_raw.size());

	int if_jpg_index = search_vct(img_raw, "ffd8");

	if (if_jpg_index > -1) {
		printf("This is JPG file.\n");
	}
	else {
		printf("This is not JPG file.\n");
		return 1;
	}

	int end_index = search_vct(img_raw, "ffd9");

	VECTARR(BYTE) ffdb = get_segments(img_raw, "ffdb");
	VECTARR(BYTE) ffc4 = get_segments(img_raw, "ffc4");
	VECTARR(BYTE) ffcx_sof = get_segments(img_raw, "ffc0");
	if (ffcx_sof.size() == 0) ffcx_sof = get_segments(img_raw, "ffc2");
	VECTARR(BYTE) ffda = get_segments(img_raw, "ffda");

	printf("***DHT***\n");
	vector<dht_t> dht_tables = ana_dht(ffc4);
	for (int i = 0; i < dht_tables.size(); ++i) {
		printf("%d-%d\n", dht_tables[i].tc, dht_tables[i].th);
		for (int j = 0; j < dht_tables[i].haffs.size(); ++j) {
			haff_t t = dht_tables[i].haffs[j];
			printf("haff:%s, rl:%d, val:%d\n", t.haff.c_str(), t.r_len, t.val);
		}
	}

	printf("***SOF***\n");
	for (int i = 0; i < ffcx_sof.size(); ++i) {
		indicate(ffcx_sof[i]);
	}
	vector<sof_t> sof_tables = ana_sof(ffcx_sof);
	for (int i = 0; i < sof_tables.size(); ++i) {
		sof_t t = sof_tables[i];
		printf("Width: %dpx Height: %dpx \n", t.wx, t.hy);
		for (int j = 0; j < t.p_count; ++j) {
			printf("BLOCK-%d\n", t.pors[j].id);
			printf("X-s-rate: %d Y-s-rate: %d\n", t.pors[j].vn, t.pors[j].hn);
			printf("DQ TABLE No:%d\n", t.pors[j].dqn);
		}
	}

	printf("***SOS***\n");
	for (int i = 0; i < ffda.size(); ++i) {
		indicate(ffda[i]);
	}



	return 0;

}
