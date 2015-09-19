#ifndef JSONPARSER_H
#define JSONPARSER_H

#include <string>
#include <vector>
#include <map>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

class JSONParser {

 public:
  JSONParser();
  ~JSONParser();
  
  void parseJSONFile( const std::string & file_name );
  bool isAGoodLumi ( int run_number, int lumi );
  void printGoodLumis();

 private:

  typedef std::map < int , std::vector<std::pair < int, int > > > GoodLumiMap;
  typedef std::map < int , std::vector<std::pair < int, int > > >::iterator GoodLumiMapIterator;
  typedef std::pair < int, std::vector < std::pair < int , int > > > GoodLumiMapEntry;

  GoodLumiMap m_good_lumi_map ;
  std::string m_file_name;

  void addToMap ( int run_number , const std::vector<std::pair<int,int> > & lumi_ranges ) ;

};

#endif
