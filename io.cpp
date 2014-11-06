#include <io.h>
#include <system.h>
#include <atom.h>
#include <unitconverter.h>
#include <cstdlib>
#include "cpelapsedtimer.h"

using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::open(char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system, bool withGhosts)
{
    CPElapsedTimer::getInstance().disk().start();
    int count = system->atomCount();

    if (withGhosts) {
        for (auto &block : system->m_ghostBlocks) {
            count += block.counter;
        }
    }

    file << count << endl;
    file << "The is an optional comment line that can be empty." << endl;
    system->for_each([&](AtomBlock &block, int i) {
        float x = block.position.x[i];
        float y = block.position.y[i];
        float z = block.position.z[i];
        file << "Ar " << UnitConverter::lengthToAngstroms(x) << " " << UnitConverter::lengthToAngstroms(y) << " " << UnitConverter::lengthToAngstroms(z) << endl;
    });

    if (!withGhosts) return;

    for (auto &block : system->m_ghostBlocks) {
        for (int i = 0; i < block.counter; i++) {
            float x = block.position.x[i];
            float y = block.position.y[i];
            float z = block.position.z[i];
            file << "Ar " << UnitConverter::lengthToAngstroms(x) << " " << UnitConverter::lengthToAngstroms(y) << " " << UnitConverter::lengthToAngstroms(z) << endl;
        }
    }

    CPElapsedTimer::getInstance().disk().stop();
}
