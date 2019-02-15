#include "BRDFReader.h"
#include <experimental/filesystem>


namespace ChefDevr {

    void BRDFReader::extract_brdfFilePaths(const char *fileDirectory) {
        using namespace std::experimental::filesystem;

        if (!is_directory(fileDirectory)) {
            throw BRDFReaderError{"The directory " + std::string{fileDirectory} + " does not exist"};
        }

        for (const path &filePath : directory_iterator(fileDirectory)) {
            brdf_filePaths.push_back(filePath);
            const auto filename = filePath.filename();
            brdf_filenames.push_back(filename.string());
        }
    }

}