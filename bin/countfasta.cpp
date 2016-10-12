// convert a fastq file to fasta
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
using namespace std;

int main (int argc, char* argv[]) {

  if ( argc < 2) {
    cerr << "Usage: " << argv[0] << " <fasta file>\n";
    return 1;
  }


  ifstream fasta;
  fasta.open(argv[1]);
  if (!fasta) return 1;

  string line;

  int minln = 10000000;
  string minid;
  int maxln = 0;
  string maxid;

  int num=0;
  long totallen=0;
  string seq = "";
  string id;
  string gt = ">";

  while (getline(fasta, line)) {
    if (line.length() == 0) {
      continue;
    }

    if ((line.at(0)) == '>') {
      if (seq.length() > 0) {
        totallen += seq.length();
        num++;
        if (seq.length() > maxln) {
          maxid = id;
          maxln = seq.length();
        }
        if (seq.length() < minln) {
          minid = id;
          minln = seq.length();
        }
      }
      id = line;
      seq = "";
    }
    else {
      seq += line;
    }
  }

  if (seq.length() > 0) {
    totallen += seq.length();
    num++;
    if (seq.length() > maxln) {
      maxid = id;
      maxln = seq.length();
    }
    if (seq.length() < minln) {
      minid = id;
      minln = seq.length();
    }
  }


  std::cout << "Length: " << totallen << "\n";
  std::cout << "Number of sequences: " << num << "\n";
  std::cout << "Longest: " << maxln << " (" << maxid << ")\n";
  std::cout << "Shortest: " << minln << " (" << minid << ")\n";


  return 0;
}
