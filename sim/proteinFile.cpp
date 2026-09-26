#include "proteinFile.hpp"

std::filesystem::path getUserHomeDir() {
    #if defined(_WIN32)
    const char* home = std::getenv("USERPROFILE");
    #else
    const char* home = std::getenv("HOME");
    #endif

    if (!home) {
        return "";
    }

    return std::filesystem::path(home);
}

std::vector<std::string> listSavedProteins() {
    std::vector<std::string> prots;
    for (const auto& f : std::filesystem::directory_iterator(getUserHomeDir() / "Simulator")) {
        if (f.path().extension() == ".protein") {
            prots.push_back(f.path().stem().string());
        }
    }

    return prots;
}

std::string convertToJoinedList(std::vector<std::string> proteins) {
    std::string prots = "";
    int i = 0;
    for (auto& n : proteins) {
        if (i > 0) prots += ";";
        prots += n;
        i++;
    }

    return prots;
}