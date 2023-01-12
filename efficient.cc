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
using namespace std;

extern int errno;

std::vector<std::vector<int>> alpha = {
  {0,110,48,94},
  {110,0,118,48},
  {48,118,0,110},
  {94,48,110,0}
};

std::unordered_map<char, int> acidMap = {
    {'A',0},{'C',1},{'G',2},{'T',3}
};

int delta = 30;

int getNums(std::vector<std::string>& lines, int start);

std::vector<int> generatePVector(std::vector<std::string>& input, int start, int end);

std::string DNAProduce(std::string & base, std::vector<int>& pS);

long getTotalMemory();

//initialize a table where the dynamic programming will be processed
std::vector<std::vector<int>> buildTable(std::string& a, std::string& b);

//initialize a table to find the alignment
std::vector<std::string> findAlignment(std::string& str1, std::string& str2);

vector<vector<int>> efficientTwoVectors(string& a, string& b);

vector<vector<int>> efficientTwoVectorsHelper(vector<vector<int>> &twoVectors, int posXH, string& a, string& b);

int getTotalCost(std::vector<int>& v1, std::vector<int>& v2);

int getIndex(std::vector<int>& v1, std::vector<int>& v2);

std::vector<std::string> findAlignmentS(std::string& a, std::string& b);

std::string getString(std::vector<char>& characters);

int main(int argc, char* argv[]) {
    
  if (argc!= 3)
  {
    std::cout << "Invalid arguments!";
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

  std::string XL = DNA1.substr(0,DNA1.length()/2);
  std::string DNA1R = DNA1.substr(DNA1.length()/2);

  std::ofstream opfile;
  opfile.open(argv[2]);

  // reverse the DNA1R and note it as XR
  string XR = DNA1R;
  reverse(XR.begin(),XR.end());

  //reverse the DNA2;
  std::string DNA2REV = DNA2;
  reverse(DNA2REV.begin(), DNA2REV.end());
  // std::vector<std::vector<int>> tableR(XR.size()+1, std::vector<int>(DNA2REV.size(),0));
  // tableR = buildTable(DNA2REV, XR);

  std::vector<std::vector<int>> twoVectors(2, std::vector<int>(DNA2.length()+1,0));
  twoVectors = efficientTwoVectors(DNA2, XL);

  std::vector<int> check1 = efficientTwoVectors(DNA2, XL)[1];
  std::vector<int> check2 = efficientTwoVectors(DNA2REV, XR)[1];

  int theCost = getTotalCost(check1,check2);
  //int theIndex = getIndex(check1,check2);
  opfile << "-------------------------------------"<<std::endl;
  opfile << theCost << std::endl;
  opfile << "-------------------------------------"<<std::endl;

  std::vector<std::string> newTableString;
  newTableString = findAlignmentS(DNA1,DNA2);
  opfile << newTableString.at(0) << '\n'<< newTableString.at(1) << '\n';
    // opfile << newTableString.at(0) << std::endl << newTableString.at(1) << std::endl;

  // Check the sequence alignment we get
  double sum = 0;
  for(int i = 0; i < int(newTableString[0].size()); i++){
    if(newTableString[0][i] == '_' || newTableString[1][i] == '_'){
      sum += 30;  
    }else{
      sum += alpha[acidMap[newTableString[0][i]]][acidMap[newTableString[1][i]]];
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

std::vector<std::vector<int>> buildTable(std::string& a, std::string& b){
  int a_len = a.length();
  int b_len = b.length();
  std::vector<std::vector<int>> dp(b_len+1, std::vector<int>(a_len+1,0));
  for(int j=0; j<=b_len; j++){
    for(int i=0; i<=a_len; i++){
      if(i == 0){
        dp[j][i] = j * delta;
        continue;
      }
      if(j == 0){
        dp[j][i] = i * delta;
        continue;
      }
      dp[j][i] = std::min(std::min(dp[j-1][i] + delta, dp[j][i-1] + delta),
                          dp[j-1][i-1] + alpha[acidMap[a[i-1]]][acidMap[b[j-1]]]);
    }
  }
  return dp;
}

vector<vector<int>> efficientTwoVectors(string& a, string& b){
  vector<vector<int>> twoVectors(2,vector<int>(a.length()+1,0));
  twoVectors = efficientTwoVectorsHelper(twoVectors, 1, a, b);
  return twoVectors;
}

vector<vector<int>> efficientTwoVectorsHelper(vector<vector<int>> &twoVectors, int posXH, string& a, string& b){
  int a_len = a.length();
  int b_len = b.length();
  if(posXH == 1){
    for(int i=0; i<=a_len; i++){
      twoVectors[0][i] = i * delta;
    }
    twoVectors[1][0] = 1 * delta;
    for(int i=1; i<=a_len; i++){
      twoVectors[1][i] = std::min(std::min(twoVectors[0][i]+delta, twoVectors[1][i-1]+delta),
                                  twoVectors[0][i-1]+alpha[acidMap[b[0]]][acidMap[a[i-1]]]);
    }
  }
  if(posXH == b_len){
      return twoVectors;
  }
  //else, the posXH must between 1 and b_len
  posXH++;
  twoVectors[0] = twoVectors[1];
  twoVectors[1][0] = twoVectors[1][0] + delta;
  for(int i=1; i<=a_len; i++){
    twoVectors[1][i] = std::min(std::min(twoVectors[0][i]+delta, twoVectors[1][i-1]+delta),
                                twoVectors[0][i-1]+alpha[acidMap[b[posXH-1]]][acidMap[a[i-1]]]);
  }

  return efficientTwoVectorsHelper(twoVectors, posXH, a, b);

}

// find the last colum using iteration
std::vector<int> buildTableS(std::string& a, std::string& b){
  int a_len = a.length();
  int b_len = b.length();
  std::vector<std::vector<int>> dp(2,std::vector<int>(a_len+1,0));
  int prev = 0, cur = 1;
  // initialize
  for(int i = 0; i <= a_len; i++){
    dp[prev][i] = i * delta;
  }
  for(int j = 0; j <= b_len; j++){
    dp[cur][0] = (j+1) * delta;
    for(int i = 1; i <= a_len; i++){
    dp[cur][i] = std::min(std::min(dp[prev][i] + delta, dp[cur][i-1] + delta),
        dp[prev][i-1] + alpha[acidMap[b[j-1]]][acidMap[a[i-1]]]);
    }
    int temp = cur;
    cur = prev;
    prev = temp;
  }
  return dp[prev];
}

// 输入两个一维数组，算出最佳匹配的的最小cost，仅仅只能得到cost，前提是输入的两个一维数组是正确的
int getTotalCost(std::vector<int>& v1, std::vector<int>& v2){
  int steps = v1.size()-1;
  int cost = v1.at(0) + v2.at(steps);
  for(int i=0; i<=steps; i++){
    if(v1.at(i) + v2.at(steps-i) < cost){
      cost = v1.at(i)+v2.at(steps-i);
    }
  }
  return cost;
}

// 输入两个一维数组，得到分割子问题所对应的index，分割点（index）是输入的v1所对应的字符串，前提是输入的两个一维数组是正确的
int getIndex(std::vector<int>& v1, std::vector<int>& v2){

  int index = 0;
  int steps = v1.size()-1;
  int cost = v1.at(0) + v2.at(steps);

  for(int i=0; i<=steps; i++){
    if(v1.at(i) + v2.at(steps-i) < cost){
      cost = v1.at(i)+v2.at(steps-i);
      index = i;
    }
  }
  return index;
}

// 根据两条字符串a,b，及函数内部构建table，无需传入，以此生成字符串的数组ss，a对应ss[0], b对应ss[1]
vector<string> findAlignment(std::string& a, std::string& b){

  vector<vector<int>> opt(b.length()+1,vector<int>(a.length()+1,0));
  opt = buildTable(a,b);
  string align1, align2;
  int left, below, diag;
  int col = a.size();
  int row = b.size();

  while(row != 0 && col !=0){
    below = opt[row-1][col] + delta;
    left = opt[row][col-1] + delta;
    diag = opt[row-1][col-1] + alpha[acidMap[a[col-1]]][acidMap[b[row-1]]];
    // the box we are interested whose optimal value comes from diag
    if(diag<=left && diag<=below){
      align1 += a[col-1];
      align2 += b[row-1];
      row --;
      col --;
      continue;
    }
    // the box we are interested whose optimal value comes from left
    if(left<=below && left<=diag){
      align1 += a[col-1];
      align2 += '_';
      col --;
      continue;
    }
    // the box we are interested whose optimal value comes from below
    if(below<=left && below<=diag ){
      align1 += '_';
      align2 += b[row-1];
      row --;
      continue;
    }
  }

  if(col == 0 && row !=0){
    for(int i=row; i>0; i--){
      align1 += '_';
      align2 += b [i-1];
    }
  }

  if(row == 0 && col !=0){
    for(int i=col; i>0; i--){
      align1 += a[i-1];
      align2 += '_';
    }
  }
  reverse(align1.begin(),align1.end());
  reverse(align2.begin(),align2.end());
  std::vector<std::string> ss(2);
  ss.at(0) = align1;
  ss.at(1) = align2;
  return ss;
}

// 内存优化算法，去找到best alignment，储存在string的vector中
std::vector<std::string> findAlignmentS(std::string& a, std::string& b){

  std::vector<std::string> ss(2);

  if(a.length() <= 2 ){
    return findAlignment(a,b);
  }

  std::string aL = a.substr(0,a.length()/2);
  std::string aR = a.substr(a.length()/2);

  string aRReversed = aR;
  reverse(aRReversed.begin(), aRReversed.end());
  string bReversed = b;
  reverse(bReversed.begin(),bReversed.end());

  vector<vector<int>> twoVectorsL;
  vector<vector<int>> twoVectorsR;
  twoVectorsL =  efficientTwoVectors(b,aL);
  twoVectorsR =  efficientTwoVectors(bReversed,aRReversed);

  int index = getIndex(twoVectorsL[1],twoVectorsR[1]);

  string bL = b.substr(0,index);
  string bR = b.substr(index);

  vector<string> ssL = findAlignmentS(aL,bL);
  vector<string> ssR = findAlignmentS(aR,bR);


  ss[1] = ssL[1] + ssR[1];// for b combination
  ss[0] = ssL[0] + ssR[0];// for a combination

  return ss;

}
