#include <iostream>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <sys/resource.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>

extern int errno;

long getTotalMemory();
int getNums(std::vector<std::string>& lines, int start);
std::vector<int> generatePVector(std::vector<std::string>& input, int start, int end);
std::string DNAProduce(std::string& base, std::vector<int>& pS);
std::vector<std::vector<int>> BuildTable(std::string& a, std::string& b, std::vector<std::vector<int>>& alpha, 
                          std::unordered_map<char, int>& map, int delta);
std::vector<std::vector<char>> FindAlignment(std::string& str1, std::string& str2, std::vector<std::vector<int>>& opt, 
                                  std::vector<std::vector<int>>& alpha, std::unordered_map<char, int>& map, int delta);                                                            

int main(int argc, char* argv[]) 
{
  if (argc!= 3)
  {
    std::cout << "number of args error"<<std::endl;
    return 0;
  }
  

  std::ifstream ifs(argv[1],std::ios::in);
  if(!ifs.is_open()){
    std::cerr << "Open File Failed" << std::endl;
    return 0;
  }

  // compress the whole file into one string vector
  std::vector<std::string> items;
  std::string oneLine;

  while (getline(ifs,oneLine)) {
    items.push_back(oneLine);
  }


  // items to be processed out
  std::string S1;
  std::string S2;
  std::vector<int> pS1;
  std::vector<int> pS2;

  S1 = items.at(0);
  int posNums1 = getNums(items,1);
  pS1 = generatePVector(items,1,posNums1+1);
  
  S2 = items.at(posNums1+1);
  pS2 = generatePVector(items,posNums1+2,items.size());

  std::string DNA1 = DNAProduce(S1,pS1);
  std::string DNA2 = DNAProduce(S2,pS2);

  //Memory and time caculation
  struct timeval begin, end;
  gettimeofday(&begin, 0);

  //Mapping A,C,G,T to 0,1,2,3
  std::unordered_map<char, int> map;
  map['A'] = 0; map['C'] = 1; map['G'] = 2; map['T'] = 3; 

  //Initialize the alpha graph
  std::vector<std::vector<int>> alpha(4, std::vector<int>(4, 0));
  alpha[0][0] = 0;   alpha[0][1] = 110; alpha[0][2] = 48;  alpha[0][3] = 94;
  alpha[1][0] = 110; alpha[1][1] = 0;   alpha[1][2] = 118; alpha[1][3] = 48;
  alpha[2][0] = 48;  alpha[2][1] = 118; alpha[2][2] = 0;   alpha[2][3] = 110;
  alpha[3][0] = 94;  alpha[3][1] = 48;  alpha[3][2] = 110; alpha[3][3] = 0;

  //Gap penalty
  int delta = 30;

  //Tabulation
  std::vector<std::vector<int>> opt = BuildTable(DNA1, DNA2, alpha, map, delta);

  //Sequence alignment
  std::vector<std::vector<char>> alignment = FindAlignment(DNA1, DNA2, opt, alpha, map, delta);
  std::reverse(alignment[0].begin(), alignment[0].end());
  std::reverse(alignment[1].begin(), alignment[1].end());

  std::ofstream opfile;
  opfile.open(argv[2]);

  opfile << "min_cost = " << opt[DNA2.size()][DNA1.size()] << std::endl;

  for(auto n : alignment[0]){
    opfile << n; 
  }
  opfile << std::endl;
  for(auto n : alignment[1]){
    opfile << n; 
  }
  opfile << std::endl;

  // Check the sequence alignment we get
  double sum = 0;
  for(int i = 0; i < int(alignment[0].size()); i++){
    if(alignment[0][i] == '_' || alignment[1][i] == '_'){
      sum += 30;  
    }else{
      sum += alpha[map[alignment[0][i]]][map[alignment[1][i]]];
    }
  }
  opfile << "checksum = " << sum << std::endl;

  //Memory and time caculation
  double totalmemory = getTotalMemory();
  gettimeofday(&end, 0);
  long seconds = end.tv_sec - begin.tv_sec;
  long microseconds = end.tv_usec - begin.tv_usec;
  double totaltime = seconds*1000 + microseconds*1e-3;
  opfile << totaltime << '\n';
  opfile << totalmemory << '\n';

  ifs.close();
  opfile.close();

  return 0;
} 

std::vector<int> generatePVector(std::vector<std::string>& input, int start, int end) {
  std::vector<int> output;
  int check;
  for(int j=start; j<end; j++){
    check = std::atoi(input.at(j).c_str());
    output.push_back(check);
  }
  return output;
}

int getNums(std::vector<std::string>& items, int start) {
  int count = 0;
  for (int i = start; i < int(items.size()); i++) {
    char check = items.at(i).at(0);
    if (!isdigit(check)) {
      return count;
    }
    count++;
  }
  return count;
}

std::string DNAProduce(std::string& base, std::vector<int>& pS){
  std::string temp(base.size() << pS.size(), ' ');
  temp.replace(0, base.size(), base, 0, base.size());
  for(int i = 0; i < int(pS.size()); i++){
    std::string str = temp.substr(0, (base.size() << i) );
    std::string upd = str.substr(0,pS[i]+1) + str + str.substr(pS[i]+1);
    temp.replace(0, upd.size(), upd);
  }
  return temp;
}

long getTotalMemory(){
  struct rusage usage;
  int returnCode = getrusage(RUSAGE_SELF, &usage);
  if(returnCode == 0){
    return usage.ru_maxrss;
  }else{
    printf("error %d", errno);
    return -1;
  }
}

std::vector<std::vector<int>> BuildTable(std::string& a, std::string& b, std::vector<std::vector<int>>& alpha, 
                                          std::unordered_map<char, int>& map, int delta){
  int a_len = a.length();
  int b_len = b.length();
  std::vector<std::vector<int>> d(b_len + 1, std::vector<int>(a_len + 1, 0));

  for(int i = 0; i <= a_len; i++){
    for(int j = 0; j <= b_len; j++){
      if(i == 0){
        d[j][i] = j * delta;
        continue;
      }
      if(j == 0){
        d[j][i] = i * delta;
        continue;
      }
      d[j][i] = std::min( std::min(d[j-1][i] + delta, d[j][i-1] + delta), d[j-1][i-1] + alpha[map[a[i-1]]][map[b[j-1]]]);            
    }
  }
  return d;
}

std::vector<std::vector<char>> FindAlignment(std::string& str1, std::string& str2, std::vector<std::vector<int>>& opt, 
                                  std::vector<std::vector<int>>& alpha, std::unordered_map<char, int>& map, int delta){
  std::vector<std::vector<char>> result;
  std::vector<char> align1,align2;
  int left, below, diag;     
  int row = str1.size();
  int col = str2.size();   
  
  while( row != 0 && col != 0 ){
    left = opt[col-1][row] + delta;
    below = opt[col][row-1] + delta;
    diag = opt[col-1][row-1] + alpha[map[str2[col-1]]][map[str1[row-1]]];
    if(diag <= left && diag <= below){
      align1.push_back(str1[row-1]);
      align2.push_back(str2[col-1]);
      row--;
      col--;
      continue;
    }
    if( left <= below && left <= diag){
      align1.push_back('_');
      align2.push_back(str2[col-1]);
      col--;
      continue;
    }
    if(below <= left && below <= diag){
      align1.push_back(str1[row-1]);
      align2.push_back('_');
      row--;
      continue;
    }
  }   

  if(row == 0 && col != 0){
    for(int i = col; i > 0; i--){
      align1.push_back('_');
      align2.push_back(str2[i-1]);
    }
  }
  if(row != 0 && col == 0){
    for(int i = row; i > 0; i--){
      align1.push_back(str1[i-1]);
      align2.push_back('_');
    }
  }
  result.push_back(align1);
  result.push_back(align2);
  return result;
}