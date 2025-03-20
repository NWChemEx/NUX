#include <fstream>
#include <iostream>
#include <string>
#include <vector>

bool is_single_number(const std::string& line) {
  return line.find_first_not_of("0123456789") == std::string::npos;
}

int main() {
  std::ifstream inputFile;
  inputFile.open("h2.xyz");
  if (!inputFile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  std::vector<int> coord_line_numbers;
  std::string line;
  int line_number = 0;

  while (std::getline(inputFile, line)) {
    line_number++;

    if (is_single_number(line)) {
      int num_atoms = std::stoi(line);
      std::getline(inputFile, line);
      line_number++;

      coord_line_numbers.push_back(line_number + 1);

      for (int i = 0; i < num_atoms; i++) {
        std::getline(inputFile, line);
        line_number++;
      }
    }
  }
  inputFile.clear();
  inputFile.seekg(0);

  std::string word;
  while (inputFile >> word) {
    if (word == "comment!"){
      std::cout << word << std::endl;
    }
  }
  inputFile.close();
  // std::cout << "XYZ coordinate blocks start at line numbers:\n";
  // for (int num : coord_line_numbers) {
  //   std::cout << num << std::endl;
  // }

  return 0;
}
