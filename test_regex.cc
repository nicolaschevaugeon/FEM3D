#include <string>
#include <regex>
int main(){
std::string test(" NODES" );
  std::regex exp ("NODES");
  std::cout << test << std::endl;
  if (std::regex_match (test, exp)) {
    std::cout << "Found regex " << std::endl; 
  }
  else {
    std::cout << "Not Found regex "  << std::endl; 
  }
  std::string sequence ("regular-expression");
  
  std::regex expr ("regular(.*)");

  if (std::regex_match (sequence.begin(), sequence.end(), exp)){
    std::cout << "GA" << std::endl;
  }
  {
    std::regex re("Get|GetValue");
    std::cmatch m;
    bool b0 = std::regex_search("GetValue", m, re);  // returns true, and m[0] contains "Get"
    bool b1 = std::regex_match ("GetValue", m, re);  // returns true, and m[0] contains "GetValue"
    bool b2 = std::regex_search("GetValues", m, re); // returns true, and m[0] contains "Get"
    bool b3 = std::regex_match ("GetValues", m, re); // returns false
    std::cout << b0 << b1 << b2 << b3 << std::endl;
    throw;
  }


}
