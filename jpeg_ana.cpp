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
#define BLOCK_INDEX 12
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

void indicate_block(vector<int> block){
  cout << "\n*indicate block*" << endl;
  for(int i = 0; i < block.size(); ++i){
    if(block[i] < 10 || block[i] < 0){
      cout << " " << block[i] << ", " << flush;
    }else{
      cout << block[i] << ", " << flush;
    }
    if ((i + 1) % 8 == 0){
      cout << endl;
    }
  }
  cout << endl;
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
    //printf("dht len->%d-now %d\n", len, index);
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

vector<BYTE> set_ff_to_ff00(vector<BYTE> vect){
  for(int i = 0; i < vect.size(); ++i){
    if(cmp_bin(vect[i], "ff")){
      vect.insert(vect.begin() + i + 1, 0x00);
    }
  }
  return vect;
}

vector<BYTE> insert_img_data(vector<BYTE> img_raw, vector<BYTE> new_img_data){
  int sos_index = search_vct(img_raw, "ffda");
  int sos_len = get_segment_length(img_raw, sos_index) + 2;
  int img_data_index = sos_index + sos_len + 2;
  int img_data_end_index = search_vct(img_raw, "ffd9");

  vector<BYTE> new_img(img_raw.begin(), img_raw.begin() + img_data_index);
  new_img.insert(new_img.end(), new_img_data.begin(), new_img_data.end());
  new_img.insert(new_img.end(), img_raw.begin() + img_data_end_index, img_raw.end());

  return new_img;
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

vector<BYTE> str_to_bin(string str){
  vector<BYTE> vect;
  for(int i = 0; i < str.size(); i += 8){
    string bin_str_part = str.substr(i, 8);
    dynamic_bitset<> bs(bin_str_part);
    vect.push_back((unsigned char)bs.to_ulong());
    
  }
  return vect;
}

string to_bin_str(int n) {
  bool minus = false;
  string str, minus_str;
  if(n < 0){
    minus = true;
    n *= -1;
  }
  while (n > 0) {
    str.push_back('0' + (n & 1));
    n >>= 1;
  }
  reverse(str.begin(), str.end());
  if(!minus) return str;
  boost::dynamic_bitset<> bs(str);
  to_string(bs.flip(), minus_str);
  return minus_str;
}

string search_haff(vector<haff_t> haffs, int value, int run_length = 0){
  for(int i = 0; i < haffs.size(); ++i){
    if(haffs[i].val == value && haffs[i].r_len == run_length){
      return haffs[i].haff;
    }
  }
}

char convert_bin_to_char(string binary){
  char letter;
  if(binary.size() != 8){
    return 0x03;
  }
  bitset<8> bs(binary);
  letter = (char)bs.to_ulong();
  return letter;
}

string string_to_bit(string str){
  string bin_str;
  str.push_back(0x03);
  for(int i = 0; i < str.size(); ++i){
    bitset<8> bs(str.c_str()[i]);
    bin_str += bs.to_string();
  }
  return bin_str;
}

string read_data_from_blocks(VECTARR(int) blocks){
  string hidden_bin;
  string hidden_data;
  for(int i = 2; i < blocks.size() - 1; i += 2){
    //indicate_block(blocks[i]);
    //indicate_block(blocks[i+1]);
    int num1 = blocks[i][BLOCK_INDEX];
    int num2 = blocks[i+1][BLOCK_INDEX];
    cout << num1 << " - " << num2 << endl;

    if(num1 < num2){
      hidden_bin.push_back('1');
    }else if(num1 > num2){
      hidden_bin.push_back('0');
    }else{
      cout << "Error. hidden data may be corrupted.\n" << endl;
      exit(-1);
    }
    if(hidden_bin.size() % 8 == 0){
      cout << hidden_bin << endl;
      char hidden_letter = convert_bin_to_char(hidden_bin);
      cout << hidden_letter << endl;
      if(
          hidden_letter < 0x08 ||
          (hidden_letter > 0x0b && hidden_letter < 0x20)
      ){
        break;
      }else{
        hidden_data.push_back(hidden_letter);
      }
      hidden_bin.clear();
    }
  }
  return hidden_data;
}

vector<int> make_id_box(sof_t sof){
  vector<int> id_box;
  for(int i = 0; i < sof.pors.size(); ++i){
    int por_length = sof.pors[i].vn * sof.pors[i].hn;
    id_box.insert(id_box.end(), por_length, sof.pors[i].id);
  }
  return id_box;
}

VECTARR(int) write_data(VECTARR(int) blocks, string input_data, sof_t sof){
  string bin_str = string_to_bit(input_data);
  cout << bin_str << endl;
  if(bin_str.size() * 2 > blocks.size()) {
    cout << "err.\n over size()." << endl;
    exit(-1);
  }

  vector<int> id_box = make_id_box(sof);
  int block_index = 2;
  for(int i = 0; i < bin_str.size(); ++i){
    if(id_box[block_index % id_box.size()] != 1){
      --i;
      continue;
    }
    vector<int> block1 = blocks[block_index];
    vector<int> block2 = blocks[block_index + 1];
    //indicate_block(block1);
    //indicate_block(block2);
    if( bin_str[i] == '0'){
      if(block1[BLOCK_INDEX] - 1 <= block2[BLOCK_INDEX]){
        int two_block_avr = (block1[BLOCK_INDEX] + block2[BLOCK_INDEX]) / 2;
        block1[BLOCK_INDEX] = two_block_avr + 1;
        block2[BLOCK_INDEX] = two_block_avr - 1;
      }
    }else{
      if(block1[BLOCK_INDEX] >= block2[BLOCK_INDEX] - 1){
        int two_block_avr = (block1[BLOCK_INDEX] + block2[BLOCK_INDEX]) / 2;
        block1[BLOCK_INDEX] = two_block_avr - 1;
        block2[BLOCK_INDEX] = two_block_avr + 1;
      }
    }
    //indicate_block(block1);
    //indicate_block(block2);
    blocks[block_index] = block1;
    blocks[block_index + 1] = block2;
    cout << blocks[block_index][BLOCK_INDEX] 
      << " - " << 
      blocks[block_index + 1][BLOCK_INDEX] 
      << endl;
    block_index += 2;
    if((i + 1) % 8 == 0) cout << "----------------" << endl;
  }
  return blocks;
}

VECTARR(int) ana_block(
    int             zig_zag[],
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
  //cout << bin_str << endl;
  int dc_ac_switch = 0; // 0 is dc. 1 is ac.

  vector<int> id_box = make_id_box(sof);
  
  vector<int> block;
  dht_inf table_inf;
  vector<haff_t> haffs;
  dqt_t dqt_t;
  bool eob = false;

  while (num_index < bin_len) {
    int now_block_id_index = block_index % id_box.size();
    int now_block_id = id_box[now_block_id_index];

    //cout << "now block id: " << now_block_id << endl;
    //cout << "block index: " << block_index << endl;
    //cout << "num index in block: " << num_count_in_block << endl;

    for(int i = 0; i < sos.infs.size(); ++i){
      if(sos.infs[i].id == now_block_id){
        table_inf = sos.infs[i];
        break;
      }
    }
    //cout << 
    //  "start index: " << start_index << 
    //  " now index: " << num_index << 
    //  endl;

    if(dc_ac_switch == 0) {
      string haff = bin_str.substr(start_index, num_index - start_index);
      //cout << "haff: " << haff << endl;
      for(int i = 0; i < dhts.size(); ++i){
        if(dhts[i].tc == 0 && dhts[i].th == table_inf.dct_cd){
          haffs = dhts[i].haffs;
          break;
        }
      }
      for(int i = 0; i < haffs.size(); ++i){
        if(haffs[i].haff == haff){
          //cout << "HIT.\n " <<
          //  haff << " " << haffs[i].val << endl;
          start_index = num_index;
          int data_len = haffs[i].val;
          string data_bit = bin_str.substr(start_index, data_len);
          dynamic_bitset<> raw_data_bit(data_bit);
          int data_num;
          if(data_bit[0] == '0'){
            data_num = -1 * (int)raw_data_bit.flip().to_ulong();
          } else{
            data_num = (int)raw_data_bit.to_ulong();
          }
          //cout << haff << ": " << data_num << endl;
          block.push_back(data_num);
          ++num_count_in_block;
          dc_ac_switch = 1;
          num_index += data_len;
          start_index = num_index;
          break;
        }
      }
    }else {
      string haff = bin_str.substr(start_index, num_index - start_index);
      //cout << "haff: " << haff << endl;
      for(int i = 0; i < dhts.size(); ++i){
        if(dhts[i].tc == 1 && dhts[i].th == table_inf.act_cd){
          haffs = dhts[i].haffs;
          break;
        }
      }
      int max_len = haffs.back().haff.size();
      if(haff.size() > max_len){
        cout << "haff string over. \n not found." << endl;
        exit(-1);
      }
      for(int i = 0; i < haffs.size(); ++i){
        if(haffs[i].haff == haff){
          start_index = num_index;
          int data_len = haffs[i].val;
          int zero_len = haffs[i].r_len;

          //cout << "HIT.\n" <<
          //  haff << " " << data_len << " " << zero_len << endl;
          
          if(data_len == 0 && zero_len == 0){
            //EOB
            eob = true;
          }else if(data_len == 0 && zero_len == 15){
            //ZRL
            block.insert(block.end(), 16, 0);
            num_count_in_block += 16;
          }else{
            string data_bit = bin_str.substr(start_index, data_len);
            dynamic_bitset<> raw_data_bit(data_bit);
            //cout << data_bit << endl;
            int data_num;
            if(data_bit[0] == '0'){
              data_num = -1 * (int)raw_data_bit.flip().to_ulong();
            }else{
              data_num = (int)raw_data_bit.to_ulong();
            }
            start_index += data_len;
            
            //cout << haff << ": " << data_num << endl;
            block.insert(block.end(), zero_len, 0);
            block.push_back(data_num);
            num_count_in_block += 1 + zero_len;
            num_index += data_len;
          }
          if(eob){
            block.insert(block.end(), 64 - num_count_in_block, 0);
            num_count_in_block = 64;
            eob = false;
          }
          start_index = num_index;
          break;
        }
      }
    }
    ++num_index;
    //cout << endl;

    if(num_count_in_block == 64){
      //indicate_block(block);
      vector<int> block_box(64);
      for(int i = 0; i < 64; ++i){
        block_box[zig_zag[i]] = block[i];
      }
      if(block_index > 0){
        block_box[0] += boxes[block_index - 1][0];
      }
      //indicate_block(block_box);
      boxes.push_back(block_box);
      dc_ac_switch = 0;
      num_count_in_block = 0;
      block.clear();
      ++block_index;
    }
  }

  return boxes;
}

string block_to_img(
    int           zig_zag[],
    vector<dqt_t> dqts,
    vector<dht_t> dhts,
    sos_t         sos,
    sof_t         sof,
    VECTARR(int)  blocks
  ){
  string img_str;
  vector<int> id_box;
  dht_inf table_inf;
  int dc_ac_switch = 0; // 0 is dc, 1 is ac.

  for(int i = 0; i < sof.pors.size(); ++i){
    int por_length = sof.pors[i].vn * sof.pors[i].hn;
    id_box.insert(id_box.end(), por_length, sof.pors[i].id);
  }
  
  for(int block_index = 0; block_index < blocks.size(); ++block_index){
    vector<int> raw_block = blocks[block_index];
    vector<int> block;
    int zrl = 0;
    int zero_start_index = 0;
    bool eob = false;

    for(int i = 0; i < 64; ++i) block.push_back(raw_block[zig_zag[i]]);
    if (block_index > 0) block[0] -= blocks[block_index - 1][0];
    //indicate_block(block);

    int now_block_id = id_box[block_index % id_box.size()];
    for( int i = 0; i < sos.infs.size(); ++i){
      if(sos.infs[i].id == now_block_id){
        table_inf = sos.infs[i];
        break;
      }
    }

    vector<haff_t> dc_haffs, ac_haffs;
    for(int i = 0; i < dhts.size(); ++i){
      if(dhts[i].tc == 0 && dhts[i].th == table_inf.dct_cd){
        dc_haffs = dhts[i].haffs;
      }
      if(dhts[i].tc == 1 && dhts[i].th == table_inf.act_cd){
        ac_haffs = dhts[i].haffs;
      }
    }
    
    for(int num_index = 0; num_index < block.size(); ++num_index){

      if(num_index == 0) {
        dc_ac_switch = 0;
        string bin_num;
        int bin_num_len;

        if(block[num_index] == 0){
          bin_num_len = 0;
        }else{
          bin_num = to_bin_str(block[num_index]);
          bin_num_len = bin_num.size();
        }
        string haff = search_haff(dc_haffs, bin_num_len);
        //cout << "haff-bin_num-len" << endl;
        //cout << haff << "-" << bin_num << "-" << bin_num_len << endl;
        img_str += haff;
        if(block[num_index] != 0) img_str += bin_num;
      }
      else {
        dc_ac_switch = 1;
        if(block[num_index] == 0){
          if (zrl == 0) zero_start_index = num_index;
          ++zrl;
          eob = true;
        }else{
          if(zrl >= 16){
            string zrl_haff = search_haff(ac_haffs, 0, 15);
            while(zrl >= 16){
              //cout << "zrl: " << zrl_haff << endl;
              img_str += zrl_haff;
              zrl -= 16;
            }
            string bin_num = to_bin_str(block[num_index]);
            //cout << "zero: " << zrl << endl;
            string haff = search_haff(ac_haffs, bin_num.size(), zrl);
            img_str += haff + bin_num;
            //cout << "haff-bin_num" << endl;
            //cout << haff << "-" << bin_num << endl;
          }else{
            string bin_num = to_bin_str(block[num_index]);
            //cout << "zero: " << zrl << endl;
            string haff = search_haff(ac_haffs, bin_num.size(), zrl);
            img_str += haff + bin_num;
            //cout << "haff-bin_num" << endl;
            //cout << haff << "-" << bin_num << endl;
          }
          zrl = 0;
          eob = false;
        }
      }
    }
    if(eob){
      string eob_haff = search_haff(ac_haffs, 0, 0);
      //cout << "eob: " << eob_haff << endl;
      img_str += eob_haff;
      eob = false;
      zrl = 0;
    }
  }
  int len = img_str.size();
  if(len % 8 != 0){
    for(int i = 0; i < (8 + 1) - (i % 8); ++i) img_str.push_back('1');
  }

  return img_str;
}


int main(int argc, char* argv[]) {
  const char *filepath = argv[1];
  const char *mode = argv[2];
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
  VECTARR(int) blocks = ana_block(zig_zag, img_data, dqt_tables, dht_tables, sos_tables[0], sof_tables[0]);

  string result;
  if(strcmp(mode, "-w") == 0){
    cout << "writing mode." << endl;
    blocks = write_data(blocks, string("@nei_phi213km"), sof_tables[0]);
    cout << "written" << endl;
  }else if(strcmp(mode, "-r") == 0){
    cout << "reading mode." << endl;
    result = read_data_from_blocks(blocks);
    cout << result << endl;
    exit(0);
  }else{
    cout << "only check header." << endl;
    exit(0);
  }

  string back_img_str = block_to_img(zig_zag, dqt_tables, dht_tables, sos_tables[0], sof_tables[0], blocks);

  vector<BYTE> new_img_data = str_to_bin(back_img_str);

  new_img_data = set_ff_to_ff00(new_img_data);

  vector<BYTE> new_img_content = insert_img_data(img_raw, new_img_data);

  ofstream ofs;
  ofs.open("output.jpg", ios::out|ios::binary|ios::trunc);
  ofs.write((const char*)&new_img_content[0], new_img_content.size());
/*
  cout << cat_bin_str(img_data) << endl;
  cout << back_img_str << endl;
  cout << back_img_str.size() << endl;

  cout << (back_img_str == cat_bin_str(img_data)) << endl;
  
  indicate(str_to_bin(back_img_str));
  cout << str_to_bin(back_img_str).size() << endl;
  cout << endl;
  indicate(img_data);

  cout << (img_data == str_to_bin(back_img_str)) << endl;

*/
  return 0;
  
}
