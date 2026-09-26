#ifndef SIM_PROTEIN_FILE_HPP
#define SIM_PROTEIN_FILE_HPP

#include <filesystem>
#include <string>
#include <vector>

std::filesystem::path getUserHomeDir();
std::vector<std::string> listSavedProteins();
std::string convertToJoinedList(std::vector<std::string>);

#endif