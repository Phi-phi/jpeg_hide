#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <vector>
#include <cstring>
#include <climits>
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
  int dct_cd;;
} dht_inf;

typedef struct sos {
  int dht_inf_count;
  int s_start;
  int s_end;
  int ah;
  int al;
  vector<dht_inf> infs;
} sos_t;

typedef struct image_infomation {
  int all_size;
  int img_type;
} image_inf;

typedef struct segment_tables {
  vector<BYTE> dqt;
  vector<BYTE> dht;
  vector<BYTE> sof;
  vector<BYTE> sos;
  vector<BYTE> img_data;
} segments;

// for debug
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

vector<dht_inf> ana_sos_block(vector<BYTE> vect, int count) {
  vector<dht_inf> tables;
  int index = 2;
  for (int i = 0; i < count; ++i) {
    dht_inf t;
    t.id = vect[++index];
    int t_data = vect[++index];
    t.dct_cd = t_data >> 4;
    t.act_cd = t_data & 0x0f;
    tables.push_back(t);
  }
  return tables;
}

vector<sos_t> ana_sos(VECTARR(BYTE) vect) {
  vector<BYTE> seg;
  vector<sos_t> tables;
  for (int i = 0; i < vect.size(); ++i) {
    seg = vect[0];
    sos_t t;
    int len = (seg[0] << 8) + seg[1] - 2;
    t.dht_inf_count = seg[2];
    int after_dht_block_index = t.dht_inf_count * 2 + 3;
    t.s_start = seg[after_dht_block_index];
    t.s_end = seg[++after_dht_block_index];
    t.ah = seg[++after_dht_block_index] >> 4;
    t.al = seg[after_dht_block_index] & 0x0f;
    t.infs = ana_sos_block(seg, t.dht_inf_count);
    tables.push_back(t);
  }
  return tables;
}

int get_segment_length(vector<BYTE> vect, int seg_index) {
  return (vect[seg_index + 2] << 8) + vect[seg_index + 3] - 2;
}

vector<BYTE> get_image_data(vector<BYTE> vect, int img_mode) {
  vector<BYTE> result;
  int sos_index = search_vct(vect, "ffda");
  int sos_len = get_segment_length(vect, sos_index) + 2; // 2 is marker's lengt
  if (img_mode == 1) {
    int end_index = search_vct(vect, "ffd9");
    for (int i = sos_index + sos_len + 2; i < end_index; ++i){
      result.push_back(vect[i]);
    }
  }
  return result;
}

vector<BYTE> set_ff00_to_ff(vector<BYTE> vect) {
  for (int i = 0; i < vect.size() - 1; ++i) {
    if (cmp_bin(vect[i], "ff") && cmp_bin(vect[i + 1], "00")) {
      vect.erase(vect.begin() + i + 1, vect.begin() + i + 2);
    }
  }
  return vect;
}

string cat_bin_str(vector<BYTE> vect) {
  int vect_len = vect.size();
  string bin_str;
  for (int i = 0; i < vect_len; ++i) {
    dynamic_bitset<> bit_value(8, vect[i]);
    string bin_str_part;
    to_string(bit_value, bin_str_part);
    bin_str += bin_str_part;
  }

  return bin_str;
}

VECTARR(int) ana_block(
    vector<BYTE>    img_data,
    vector<dqt_t>   dqts,
    vector<dht_t>   dhts,
    sos_t           sos,
    sof_t           sof
  ) {
  VECTARR(int) boxes;
  string bin_str = cat_bin_str(img_data);
  int start_index = 0;
  int num_index = 0;
  int block_index = 0;
  int num_count_in_block = 0;
  int bin_len = bin_str.size();
  cout << bin_str << endl;
  int dc_ac_switch = 0; // 0 is dc. 1 is ac.

  vector<int> id_box;
  vector<int> dqt_box;
  
  vector<int> block;
  dht_inf table_inf;
  vector<haff_t> haffs;
  dqt_t dqt_t;

  for(int i = 0; i < sof.pors.size(); ++i){
    int por_length = sof.pors[i].vn * sof.pors[i].hn;
    id_box.insert(id_box.end(), por_length, sof.pors[i].id);
    dqt_box.push_back(sof.pors[i].dqn);
  }

  while (num_index < bin_len) {
    int now_block_id_index = block_index % id_box.size();
    int now_block_id = id_box[now_block_id_index];

    for(int i = 0; i < sos.infs.size(); ++i){
      if(sos.infs[i].id == now_block_id){
        table_inf = sos.infs[i];
        break;
      }
    }

    if(dc_ac_switch == 0) {
      printf("DC\n");
      string haff = bin_str.substr(start_index, num_index);
      for(int i = 0; i < dhts.size(); ++i){
        if(dhts[i].tc == 0 && dhts[i].tc == now_block_id){
          haffs = dhts[i].haffs;
          break;
        }
      }
      for(int i = 0; i < haffs.size(); ++i){
        if(haffs[i].haff == haff){
          start_index += num_index;
          int data_len = haffs[i].val;
          string data_bit = bin_str.substr(start_index, data_len);
          dynamic_bitset<> raw_data_bit(data_bit);
          int data_num;
          if(data_bit[0] == 1){
            data_num = -1 * (int)raw_data_bit.flip().to_ulong();
          } else{
            data_num = (int)raw_data_bit.to_ulong();
          }
          cout << haff << ": " << data_num << endl;
          block.push_back(data_num);
          ++num_count_in_block;
          dc_ac_switch = 1;
          start_index += data_len;
          break;
        }
        ++num_index;
      }
    }else {
      printf("AC\n");
      string haff = bin_str.substr(start_index, num_index);
      for(int i = 0; i < dhts.size(); ++i){
        if(dhts[i].tc == 1 && dhts[i].tc == now_block_id){
          haffs = dhts[i].haffs;
          break;
        }
      }
      for(int i = 0; i < haffs.size(); ++i){
        if(haffs[i].haff == haff){
          start_index += num_index;
          int data_len = haffs[i].val;
          int zero_len = haffs[i].r_len;

        }
      }
    }
  }

  return boxes;
}

int main(int argc, char* argv[]) {
  const char *filepath = argv[1];
  vector<BYTE> img_raw;
  image_inf inf;
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

  inf.all_size = img_raw.size();

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
  if (ffcx_sof.size() == 0) {
    inf.img_type = 2;
    ffcx_sof = get_segments(img_raw, "ffc2");
  }

  VECTARR(BYTE) ffda = get_segments(img_raw, "ffda");

  vector<dqt_t> dqt_tables = ana_dqt(ffdb);

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
  sof_t t = sof_tables[0];
  printf("Width: %dpx Height: %dpx \n", t.wx, t.hy);
  for (int j = 0; j < t.p_count; ++j) {
    printf("BLOCK-%d\n", t.pors[j].id);
    printf("X-s-rate: %d Y-s-rate: %d\n", t.pors[j].vn, t.pors[j].hn);
    printf("DQ TABLE No:%d\n", t.pors[j].dqn);
  }

  printf("***SOS***\n");
  for (int i = 0; i < ffda.size(); ++i) {
    indicate(ffda[i]);
  }
  vector<sos_t> sos_tables = ana_sos(ffda);
  for (int i = 0; i < sos_tables.size(); ++i) {
    sos_t t = sos_tables[i];
    for(int j = 0; j < t.dht_inf_count; ++j) {
      dht_inf inf = t.infs[j];
      printf("id: %d, ac: %d, dc: %d\n", inf.id, inf.act_cd, inf.dct_cd);
    }
    printf("s_start: %d, s_end: %d\n", t.s_start, t.s_end);
    printf("Before scan: %d\n", t.ah);
    printf("Now scan: %d\n", t.al);
  }

  printf("***IMAGE DATA***\n");
  vector<BYTE> img_data = get_image_data(img_raw, 1);
  printf("len- %d\n", (int)img_data.size());

  img_data = set_ff00_to_ff(img_data);
  printf("Img data length %d\n", (int)img_data.size());
  ana_block(img_data, dqt_tables, dht_tables, sos_tables, sof_tables[0]);

  return 0;

}
